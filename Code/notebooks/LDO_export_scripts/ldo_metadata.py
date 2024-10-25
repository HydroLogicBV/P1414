# -*- coding: utf-8 -*-
"""
Generates the file required for the LDO database based on the D-HYDRO output files.

Created on Fri December 1th 2023
@author: Hydrologic
"""

import pandas as pd
import numpy as np
import netCDF4 as nc
import logging, re, os, json
from datetime import datetime, timedelta


class MetadataGenerator():

    def __init__(self, model_dir:str, file_locations:dict, existing_df:pd.DataFrame=None) -> None:
        # Initialize input and output directory
        self.model_dir = model_dir
        if os.path.exists(self.model_dir) == False:
            raise Exception("The specified model dir does not exist")
        
        required_locations = ['fm_dir', 'output_dir', 'mdu_file', 'map_file', 'his_file', 'structures_file', 'dia_file']
        for key in required_locations:
            if key not in file_locations.keys():
                raise Exception(f"Missing {key} in file_locations input, keys found are: {file_locations}")
        self.file_locations = file_locations

        self.model_output_dir = os.path.join(model_dir, file_locations['fm_dir'], file_locations['output_dir'])
        if os.path.exists(self.model_output_dir) == False:
            raise Exception(f"The output directory {self.model_output_dir} could not be found, are you sure your input is correct and that you ran the model?")
    
        self.write_output_dir = os.path.join(self.model_output_dir, 'LDO')
        if os.path.exists(self.write_output_dir) == False:
            os.mkdir(self.write_output_dir)


        # Initialize dictionary for metadata storage
        self.dict_meta = {}
        # Intialize dataframe for metadata storage
        if existing_df is None:
            self.df_ = pd.DataFrame()
        else:
            self.df_ = existing_df
        self.configure_logging()

    def execute(self):
        case = os.path.basename(self.model_dir)

        # Step 1: initialize dataframe
        if self.df_.empty:
            if not self.intialize_dataframe():
                raise Exception('Error when initializing the dataframe')

        # Step 2: Add a new row/scenario to the dataframe
        self.df_.loc[len(self.df_)] = np.nan

        # Step 3: Load the structures.ini file
        self.load_structures_ini(os.path.join(self.model_dir, self.file_locations['fm_dir'], self.file_locations['structures_file']))

        # Step 4: Open the NetCDF DFM_his file
        self.load_dfm_his(os.path.join(self.model_output_dir, self.file_locations['his_file']))

        # Step 5: Open the mdu dfm file
        self.load_dfm_mdu(os.path.join(self.model_dir, self.file_locations['fm_dir'], self.file_locations['mdu_file']))

        # Step 6: Open the DFM dia file
        self.load_dfm_dia(os.path.join(self.model_output_dir, self.file_locations['dia_file']))

        # Step 7: Get scenario information
        self.get_information_scenario(case)

        # Step 8: Get locatie information
        self.get_information_locatie(case)

        # Step 9: Get doorbraak information
        self.get_information_doorbraak(case)

        # Step 10: Get buitenwater information
        self.get_information_buitenwater(case)

        # Step 11: Get model information
        self.get_information_model(case)

        # Generate output file
        self.df_out = pd.DataFrame([[col[0] for col in self.df_.columns], 
                               [col[1] for col in self.df_.columns]] + self.df_.values.tolist(), columns=self.df_.columns)
        return self.df_out

    def intialize_dataframe(self):
        
        # Initialize headers to be logged regarding scenario
        self.logger.debug('Initializing header information regarding scenario')
        self.dict_meta["Scenario"] = [
                'Scenario Identificatie',
                'Scenarionaam',
                'Scenariodatum',
                'Projectnaam',
                'Eigenaar overstromingsinformatie',
                'Beschrijving scenario',
                'Versie resultaat',
                'Motivatie rekenmethode',
                'Doel',
                'Berekeningsmethode',
                'Houdbaarheid scenario' ]
        
        # Initialize headers to be logged regarding location
        self.logger.debug('Initializing header information regarding location')
        self.dict_meta["Location"] = [
                'x-coordinaten doorbraaklocatie',
                'y-coordinaten doorbraaklocatie',
                'Naam buitenwater',
                'Naam waterkering',
                'Naam doorbraaklocatie',
                'Gebiedsnaam',
                'Subgebiedsnaam' ]

        # Initialize headers to be logged regarding doorbraak
        self.logger.debug('Initializing header information regarding doorbraak')
        self.dict_meta["Bres"] = [
                'overstromingskans - Norm (dijktraject)',
                'overstromingskans - DPV referentie (dijktraject)',
                'overstromingskans - DPV referentie (trajectdeel)',
                'overstromingskans - beheerdersoordeel (trajectdeel)',
                'Keringbeheerder',
                'Materiaal kering',
                'Lengte dijkringdeel',
                'Bresdiepte',
                'Duur bresgroei in verticale richting',
                'Initiele bresbreedte',
                'Methode bresgroei',
                'Startmoment bresgroei',
                'Maximale bresbreedte',
                'Kritieke stroomsnelheid (Uc)',
                'bresgroeifactor 1 (f1)',
                'bresgroeifactor 2 (f2)',
                'Afvoer coefficient (Ce)',
                'Initial Crest [m+NAP]',
                'Rivierknoop',
                'Boezemknoop',
                'Lowest crest',
                'Gridhoogte',
                'Maximaal bresdebiet',
                'Maximale instroom',
                'Wieldiepte',
                'Maximale waterstand']

        # Initialize headers to be logged regarding buitenwater
        self.logger.debug('Initializing header information regarding buitenwater')
        self.dict_meta["Buitenwater"] = [
                'Buitenwatertype',
                'Maximale buitenwaterstand',
                'Stormvloedkering open',
                'Stormvloedkering gesloten',
                'Compartimentering van de boezem',
                'Gemiddeld meerpeil',
                'Eigenschappen getij',
                'Piekduur',
                'Stormduur',
                'Timeshift',
                'Debiet randvoorwaarden(rvw) locatie',
                'Waterstand randvoorwaarden(rvw) locatie',
                'Overschrijdingsfrequentie',
                'Maximale afvoer Rijn',
                'Maximale afvoer Maas' ]

        # Initialize headers to be logged regarding model
        self.logger.debug('Initializing header information regarding model')
        self.dict_meta["Model"] = [
            'Datum modelschematisatie',
            'Variant beschrijving',
            'Bodemhoogte model',
            'Ruwheid model',
            'Start berekening',
            'Einde berekening',
            'Rekenduur',
            'Start simulatie',
            'Einde simulatie',
            'Duur',
            'Modelversie',
            'Modelleersoftware',
            'Modelresolutie' ]

        # Initialize headers to be logged regarding overige
        self.logger.debug('Initializing header information regarding resterende')
        self.dict_meta["Resterende"] = [
            'Regionale keringen of hoge lijnelementen standzeker',
            'Overige opmerkingen'
        ]

        # Initialize headers to be logged regarding bestanden
        self.logger.debug('Initializing header information regarding bestanden')
        self.dict_meta["Bestanden"] = [
            'Maximale stroomsnelheid',
            'Animatie waterdiepte',
            'Maximale waterdiepte',
            'Rapportage',
            'Bathymetrie',
            '3Di simulatie resultaat'
        ]

        # Generate dataframe using the header information
        self.logger.debug('Generating dataframe using header information')
        cols = [(section, field)
            for section in self.dict_meta.keys()
            for field in self.dict_meta[section]
            ]
        self.df_ = pd.DataFrame(columns=cols, dtype='object')
        return True

    def write(self, write_as_excel:bool=False, output_file:str=None):
        if write_as_excel:
            if output_file == None:
                output_file = os.path.join(self.write_output_dir, 'LDO_metadata.xlsx')
            try:
                import openpyxl
            except ModuleNotFoundError:
                print("The 'openpyxl' module is not installed. Install it using 'conda install openpyxl'. If you do not want to install it, you can also set the write_as_excel parameter to False.")
            self.df_out.to_excel(output_file, index=False)
        else:                 
            if output_file == None:
                output_file = os.path.join(self.write_output_dir, 'LDO_metadata.csv')
            self.df_out.to_csv(sep=";", path_or_buf=output_file, index=False)


        self.logger.info(f"Succesfully wrote the LDO metadata to: {output_file}")
    
    # --- Get information regarding subsection --- 

    def get_information_scenario(self, case):
        self.logger.debug('Obtain scenario information for case')

        # Column A: fill in name of folder as being the scenario identification
        self.df_.iloc[-1, 0] = self.get_scenarioname(case)

        # Column C: fill in today's date as being the scenario date
        self.df_.iloc[-1, 2] = datetime.now().strftime('%Y-%m-%d')

        # Column D: fill in projectnaam
        # self.df_.iloc[-1, 3] = 'Overstromingsscenario\'s WSVV 2023'

        # Column E: fill in owner based on table
        # self.df_.iloc[-1, 4] = int(38)

        # Column G: fill in scenario version
        # self.df_.iloc[-1, 6] = float(1.0)

        # Column H: fill in motivation berekeningsmethode
        # self.df_.iloc[-1, 7] = 'Meeste detail, meenemen bresgroei, initiele waterstand meenemen, effect hoogte achterland op overstroming modelleren.'

        # Column I: fill in scenario objective
        # self.df_.iloc[-1, 8] = 'Bepalen overstroomt oppervlak als gevolg van dijkdoorbraak.'

        # Column J: fill in berekeningsmethode based on table
        # self.df_.iloc[-1, 9] = int(3)
        
        self.logger.debug('Obtain scenario information for case completed')
        return True
    
    def get_information_locatie(self, case):
        self.logger.debug('Obtain locatie information for case')

        for line in self.df_structures:
            if line.lower().startswith('StartLocationY'.lower()):
                # column M: y-coordinaten doorbraaklocatie
                key, value = map(str.strip, line.split('='))
                self.df_.iloc[-1, 12] = float(value)
                continue
               
            if line.lower().startswith('StartLocationX'.lower()):
                # column L: x-coordinaten doorbraaklocatie
                key, value = map(str.strip, line.split('='))
                self.df_.iloc[-1, 11] = float(value) 
                continue

        self.logger.debug('Obtain locatie information for case completed')
        return True

    def get_information_doorbraak(self, case):
        self.logger.debug('Obtain doorbraak information for case')

        # Column X: methode bresgroei
        self.df_.iloc[-1, 23] = int(1)

        # Columns to nvt
        # self.df_.iloc[-1, 36] = self.df_.iloc[-1, 37] = 'nvt'

        # Column AW: fill in overschrijdingsfrequentie based on foldername
        # self.df_.iloc[-1, 56] = self.get_overschrijdingsfrequentie(case)

        for line in self.df_structures:
            if line.lower().startswith('TimeToBreachToMaximumDepth'.lower()):
                # column AA: Duur bresgroei in verticale richting
                key, value = map(str.strip, line.split('='))
                self.df_.iloc[-1, 26] = self.convert_seconds_to_date(value)
                continue

            if line.lower().startswith('BreachWidthIni'.lower()):
                # column AB: Initiele bresbreedte
                key, value = map(str.strip, line.split('='))
                self.df_.iloc[-1, 27] = float(value)
                continue

            if line.lower().startswith('T0'.lower()):
                # column AD: Startmoment bresgroei
                key, value = map(str.strip, line.split('='))
                self.t0 = int(float(value))
                continue

            if line.lower().startswith('Ucrit'.lower()):
                # column AF: Kritieke stroomsnelheid (Uc)
                key, value = map(str.strip, line.split('='))
                self.df_.iloc[-1, 31] = float(value)
                continue

            if line.lower().startswith('F1'.lower()):
                # column AG: bresgroeifactor 1 (f1)
                key, value = map(str.strip, line.split('='))
                self.df_.iloc[-1, 32] = float(value)
                continue

            if line.lower().startswith('F2'.lower()):
                # column AH: bresgroeifactor 2 (f2)
                key, value = map(str.strip, line.split('='))
                self.df_.iloc[-1, 33] = float(value)
                continue
                
            if line.lower().startswith('CrestLevelIni'.lower()):
                # column AJ: Initial Crest [m+NAP]
                key, value = map(str.strip, line.split('='))
                self.df_.iloc[-1, 35] = float(value)
                continue

            if line.lower().startswith('CrestLevelMin'.lower()):
                # column AM: Lowest crest
                key, value = map(str.strip, line.split('='))
                self.df_.iloc[-1, 38] = float(value)

        # Column X: Materiaal kering
        self.df_.iloc[-1, 23] = self.get_materiaal_kering(self.df_.iloc[-1, 31])

        if 'dambreak_crest_level' in self.df_dfm_his.variables:
            # Column Z: Bresdiepte
            mask_crestlevel = self.df_dfm_his.variables['dambreak_crest_level'][:].data[:,0].flatten()
            self.df_.iloc[-1, 25] = float(max(mask_crestlevel) - min(mask_crestlevel))

        if 'dambreak_crest_width' in self.df_dfm_his.variables:
            # Column AE: Maximale bresbreedte
            mask_crestwidth = self.df_dfm_his.variables['dambreak_crest_width'][:].data[:,0].flatten()
            self.df_.iloc[-1, 30] = int(max(mask_crestwidth))

        if 'dambreak_crest_level' in self.df_dfm_his.variables:
            # column AN: Gridhoogte
            mask_crestlevel = self.df_dfm_his.variables['dambreak_crest_level'][:].data[:,0].flatten()
            self.df_.iloc[-1, 39] = float(min(mask_crestlevel))

        if 'dambreak_discharge' in self.df_dfm_his.variables:
            # column AO: Maximaal bresdebiet
            values = self.df_dfm_his.variables['dambreak_discharge'][:].data[:,0].flatten()
            mask_ddischarge = np.abs(values[values != -999.])
            self.df_.iloc[-1, 40] = float(max(mask_ddischarge))

        if 'dambreak_cumulative_discharge' in self.df_dfm_his.variables:
            # column AP: Maximale instroom
            values = self.df_dfm_his.variables['dambreak_cumulative_discharge'][:].data[:,0].flatten()
            mask_cumdischarge = np.abs(values[values != -999.])
            self.df_.iloc[-1, 41] = float(max(mask_cumdischarge))

        self.logger.debug('Obtain doorbraak information for case completed')
        return True
    
    def get_information_buitenwater(self, case):
        self.logger.debug('Obtain buitenwater information for case')

        if 'dambreak_s1up' in self.df_dfm_his.variables:
            # column AT: Maximale buitenwaterstand
            mask_dams1up = self.df_dfm_his.variables['dambreak_s1up'][:].data[:,0].flatten()
            self.df_.iloc[-1, 45] = float(max(mask_dams1up))

        self.logger.debug('Obtain buitenwater information for case completed')
        return True

    def get_information_model(self, case):
        self.logger.debug('Obtain model information for case')
        
        for line in self.dfm_mdu:
                    
            if line.lower().startswith('Tunit'.lower()):
                # Extract tunit
                tunit = line.split('=')[1].strip().split('#')[0].strip()

            if line.lower().startswith('TStart'.lower()):
                # Extract tstart
                tstart = self.convert_value_to_seconds( int(line.split('=')[1].strip().split('#')[0].strip()), tunit )

            if line.lower().startswith('TStop'.lower()):
                # Extract tstop
                tstop = self.convert_value_to_seconds( int(line.split('=')[1].strip().split('#')[0].strip()), tunit)

            if line.lower().startswith('RefDate'.lower()):
                # Extract refdate
                refdate = pd.to_datetime(line.split('=')[1].strip().split('#')[0].strip(), format="%Y%m%d")

            if line.startswith('Version'):
                # Extract version
                self.df_.iloc[-1, 69] = line.split('=')[1].strip().split('#')[0].strip()

        # Column AD: startmoment bresgroei
        self.df_.iloc[-1, 29] = self.convert_seconds_to_date(self.t0 - tstart)

        # Column BO: start simulatie
        self.df_.iloc[-1, 66] = (refdate + timedelta(seconds=tstart)).strftime('%d/%m/%Y %H:%M')

        # Column BP: stop simulatie
        self.df_.iloc[-1, 67] = (refdate + timedelta(seconds=tstop)).strftime('%d/%m/%Y %H:%M')

        # Column BQ: duur simulatie
        self.df_.iloc[-1, 68] = self.convert_seconds_to_date(((refdate + timedelta(seconds=tstop)) - (refdate + timedelta(seconds=tstart))).total_seconds())

        # Column BH: fill in date on which model has been made
        # self.df_.iloc[-1, 59] = datetime(2024).strftime('%Y-%m-%d')

        # Column BI: fill in variant beschrijving
        # self.df_.iloc[-1, 60] = 'Betreft falen van een primaire kering.'

        # Column BS: fill in variant beschrijving
        # self.df_.iloc[-1, 70] = 'D-HYDRO 1D2D 2024.03'

        comp_start = [msg for msg in self.dfm_dia if 'Computation started' in msg]
        comp_finis = [msg for msg in self.dfm_dia if 'Computation finished' in msg]
        
        # Column BL: computation started
        if comp_start:
            comp_start = datetime.strptime(comp_start[0].split(':', 2)[-1].strip(), '%H:%M:%S, %d-%m-%Y')
            self.df_.iloc[-1, 63] = comp_start.strftime('%d/%m/%Y %H:%M')
        # Column BM: computation ended
        if comp_finis:
            comp_finis = datetime.strptime(comp_finis[0].split(':', 2)[-1].strip(), '%H:%M:%S, %d-%m-%Y')
            self.df_.iloc[-1, 64] = comp_finis.strftime('%d/%m/%Y %H:%M')

        # Column BN: duur simulatie
        if comp_start and comp_finis:
            self.df_.iloc[-1, 65] = self.convert_seconds_to_date((comp_finis - comp_start).total_seconds())

        self.logger.debug('Obtain model information for case completed')
        return True

    # --- Loading functions --- 

    def load_structures_ini(self, filename):
        # Create an empty dictionary to store section data
        data, current_section = {}, None 

        # Check if structures.ini exists
        if os.path.exists(filename) == False:
            raise Exception(f"The structures file can not be found, expected it at {filename}")
        
        # Open the ini file
        with open(filename, 'r') as file:
            # Read the file line by line
            for line in file:
                line = line.strip()  # Remove leading/trailing whitespace

                # Check if the line contains a section header (e.g., [SectionName])
                if line.startswith('[') and line.endswith(']'):
                    # Extract the section name
                    current_section = line[1:-1]

                    # Check if the section doesn't exist in the data dictionary
                    if current_section not in data:
                        # Create an empty list to store lines for this section
                        data[current_section] = []
                
                # Check if the line belongs to a section (i.e., not a comment or empty line)
                elif current_section:
                    # Append the line to the list corresponding to the current section
                    data[current_section].append(line)
        self.df_structures = []

        # Check if there are one or two dambreaks
        if 'Structure' in data:                                
            for line in data['Structure'][-44:-22]:
                if line.lower().startswith('type'.lower()):
                    # column: Check if type column is dambreak
                    key, value = map(str.strip, line.split('='))
                    if value == 'dambreak':
                        # Eerste informatieblok over de bres in de primaire kering
                        self.df_structures = data['Structure'][-44:-22]
        if not self.df_structures: 
            self.df_structures = data['Structure'][-22:]

        # Return the dictionary containing section data
        return True

    def load_dfm_his(self, filename):
        if not os.path.exists(filename):
            raise Exception(f"The his file can not be found, expected it at {filename}")
        self.df_dfm_his = nc.Dataset(filename, 'r')
        return True

    
    def load_dfm_mdu(self, filename):
        if os.path.exists(filename) == False:
            raise Exception(f"The mdu file can not be found, expected it at {filename}")
        with open(filename, 'r') as file:
            self.dfm_mdu = file.readlines()
        return True
    
    def load_dfm_dia(self, filename):
        if os.path.exists(filename) == False:
            raise Exception(f"The dia file can not be found, expected it at {filename}")
        with open(filename, 'r') as file:
            self.dfm_dia = file.readlines()
        return True


    # --- Helper functions --- 
    
    def configure_logging(self):
        # Args:
        # Returns:
        #   logger: logger for the current module
        logging.basicConfig(
            level=logging.INFO,  
            # format='%(asctime)s - %(levelname)s - %(message)s', 
            handlers=[
                logging.StreamHandler() 
            ]
        )
        self.logger = logging.getLogger(__name__)

    
    def convert_value_to_seconds(self, value: int, tunit: str):
        # Args:
        #   value (int): value in the unit of measurment
        #   tunit (str): unit of measurement being one of D, H, M or S
        # Returns:
        #   value (int): value converted to seconds based on the unit of measurement
        if tunit == 'D':
            return value * 86400  # 1 day has 86400 seconds
        elif tunit == 'H':
            return value * 3600  # 1 hour has 3600 seconds
        elif tunit == 'M':
            return value * 60  # 1 minute has 60 seconds
        elif tunit == 'S':
            return value  # Value is already in seconds
        else:
            raise ValueError("Invalid time unit. Use 'D', 'H', 'M', or 'S'.")
        
    def convert_seconds_to_date(self, value: int):
        # Args:
        #   value (int): a value providing an integer value in seconds
        # Returns:
        #   date: a date in the format DD HH:MM
        days, remainder = divmod(int(float(value)), 86400)  # 86400 seconds in a day
        hours, remainder = divmod(remainder, 3600)  # 3600 seconds in an hour
        minutes, seconds = divmod(remainder, 60)  # 60 seconds in a minute
        return f"{days} d {hours:02}:{minutes:02}"

    # --- Transformation functions ---

    def get_scenarioname(self, folder: str):
        # Args:
        #   folder (str): the name of the folder of the scenario 
        # Returns:
        #   freq (int): the scenarioname of the scenario
        match = re.compile(r'^(.*?)_?T(\d+)?$').search(folder)
        if match:
            name = match.group(1).replace('_',' ')
        else:
            name = folder.replace('_',' ')
        return name

    def get_overschrijdingsfrequentie(self, folder: str):
        # Args:
        #   folder (str): the name of the folder of the scenario 
        # Returns:
        #   freq (int): the overschrijdingsfrequentie of the scenario
        match = re.compile(r'T(\d+)').search(folder)
        if match:
            freq = match.group(1)
        else:
            freq = 1000
        return freq

    def get_materiaal_kering(self, ks: int):
        # Args:
        #   ks (int): kritieke stroomsnelheid [m/s], waarbij bresgroei wordt gestart
        # Returns:
        #   mk (str): materiaal waaruit de kering is opgebouwd
        if ks == 0.2 : return 'zand'
        if ks == 0.5 : return 'klei'


if __name__ == "__main__":
    # Create an instance of MetadataGenerator    
    file_locations = {
        'fm_dir': "dflowfm",
        'output_dirname': "output",
        'mdu_file': "DFM.mdu",
        'map_file': "map.nc",
        'his_file': "DFM_his.nc",
        'structures_file': "structures.ini",
        'dia_file': "DFM.dia",
    }
    metadatagenerator = MetadataGenerator(r"C:\Work\Projects\HL-P24050_ROI\05_Analysis\SAS_workdir\run_model\test_dambreach\nothing", file_locations)

    # Call the execute method which in turn calls set_headers
    metadatagenerator.execute()

