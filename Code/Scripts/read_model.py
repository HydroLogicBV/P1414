from pathlib import Path

from hydrolib.core.io.dimr.models import DIMR

model_path = Path(r"D:\Work\Project\P1414\Models\Combined\v6_read\dimr_config.xml")
dimr_model = DIMR(filepath=model_path)
fm_component = dimr_model.component[0]
dimr_model.component[0].model.geometry.dxmin1d = 10
# print(fm_component)

dimr_model.filepath = Path(r"D:\Work\Project\P1414\Models\Combined\v6_write\dimr_config.xml")
dimr_model.save(recurse=True)
