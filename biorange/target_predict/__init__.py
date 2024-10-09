## 成分靶点的逻辑
from biorange.target_predict.mol_target.chembl_final import (
    ChEMBLTargetScraper,
    chembl_smiles_target,
)

from biorange.target_predict.mol_target.stitch_inchikey import (
    TCMDataProcessor as StichTargetIchikeyScraper,
    stich_inchikey_target,
)

# TODO 还是文件夹逻辑，不统一待修改
# from biorange.target_predict.mol_target.stitch import (
#     StichTargetSmilesScraper,
#     stich_smiles_target,
# )


from biorange.target_predict.mol_target.tcmsp_inchikey import (
    TCMSPTargetScraper as TCMSPTargetInchikeyScraper,
    tcmsp_inchikey_target,
)

from biorange.target_predict.mol_target.tcmsp_smiles import (
    TCMSPTargetScraper as TCMSPTargetSmilesScraper,
    tcmsp_smiles_target,
)

## 疾病靶点逻辑
from biorange.target_predict.disease_target.omim import (
    OmimDiseaseScraper,
    omim_disease_target,
)

from biorange.target_predict.disease_target.ttd import (
    TTDDiseaseScraper,
    ttd_disease_target,
)

from biorange.target_predict.disease_target.genecards import (
    GenecardsDiseaseScraper,
    genecards_disease_target,
)

from biorange.target_predict.data_processing.ingredient_input import admet_filter
