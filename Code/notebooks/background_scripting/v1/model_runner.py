from notebooks.background_scripting.v1.widget_styling import WidgetStyling
import os
import subprocess
import ipywidgets as ipy
from IPython.display import display, clear_output

class ModelRunner(WidgetStyling):
    def __init__(self, path):
        self.model_folder = path
        self.run_dimr_path = os.path.join(self.model_folder, 'run.bat')

        self.set_default_layout_and_styling()

        self.progess_bar = ipy.FloatProgress(
            min = 0,
            max = 100,
            value = 0,
            layout = self.progress_layout,
            style = self.progress_style
        )
        
    def display_widgets(self):
        button = ipy.Button(description="Run model!", style = self.button_style, layout = self.button_layout)
        output = ipy.Output()
        display(button, output)
        button.on_click(self.run_model)

    def read_model_progress(self):
        for index_line in reversed(range(len(self.full_logs))):
            if '%' in self.full_logs[index_line]:
                valid_percentage = self.update_model_percentage(self.full_logs[index_line])
                if valid_percentage:
                    break
        if len(self.percentage_list) > 1:
            self.initializing = False
            self.running = True

    def update_model_percentage(self, line):
        line_before_percentage = line.split('%')[0]
        percentage = line_before_percentage.split(' ')[-1]
        try:
            percentage = float(percentage)
            if percentage not in self.percentage_list:
                self.percentage_list.append(percentage)

            if len(self.percentage_list) > 1:
                self.percentage = percentage
            return True
        except:
            return False
            
    def run_model(self, b):
        nr_log_lines = 4
        self.full_logs = []
        
        self.initializing = False
        self.running = False
        self.percentage = 0.0
        self.percentage_list = []
        self.done = False

        initializing_out = ipy.HTML("")
        logging_title = ipy.HTML(f"<p>Last {nr_log_lines} logging lines: </p>", style = self.item_style, layout = self.item_layout)
        logging = ipy.HTML("", style = self.item_style, layout = self.item_layout)

        output_to_display = [self.progess_bar, initializing_out, logging_title, logging]
        items = ipy.VBox(children=output_to_display, layout=self.box_layout, style = self.box_style)
        display(items)

        with subprocess.Popen(self.run_dimr_path, cwd = self.model_folder ,  stdout=subprocess.PIPE, bufsize=1, universal_newlines=True) as p:
            self.initializing = True
            self.last_checked_line = 0
            for b in p.stdout:
                self.full_logs.append(b)
                new_logging = ""
                for line in self.full_logs[-nr_log_lines:]:
                    if "0 nr of dambreak links" in line:
                        p.kill()
                        raise Exception("Dambreak is invalid, retry setting a dambreak")
                    new_logging += f"<p>{line}</p>"
                    self.read_model_progress()
                logging.value = new_logging
                self.progess_bar.value = self.percentage
                if self.initializing:
                    initializing_out.value = "<p>Initializing model...</p>"
                else:
                    initializing_out.value = ""
        if self.percentage > 99:
            self.done = True
            logging_title.value = "<p>Model finished!</p>"
        logging.value = ""
        
        if p.returncode != 0:
            raise subprocess.CalledProcessError(p.returncode, p.args)