%% Analysis 2 subject selection
%1.find subj that don't have dose 1 and remove them
%2.separate dose 1 placs into groups of dose3-drug and dose3-plac
subj_drug_dose3 = find(~cellfun(@isempty, ERSP.Drug_dose3));
subj_plac_dose3 = find(~cellfun(@isempty, ERSP.Plac_dose3));
subj_plac_single = find(~cellfun(@isempty, ERSP.Plac_single));
subj_all_dose3 = [subj_drug_dose3 subj_plac_dose3];
subj_missing_single = setdiff(subj_all_dose3, subj_plac_single);

ERSP.Drug_dose3(subj_missing_single) = [];
ERSP.Plac_dose3(subj_missing_single) = [];

ERSP.Plac_single_d3Drug(subj_drug_dose3) = ERSP.Plac_single(subj_drug_dose3);
ERSP.Plac_single_d3Plac(subj_plac_dose3) = ERSP.Plac_single(subj_plac_dose3);
%% Analysis 2 Qmoney contrast vector creation
all_PFC_Qmoney_gamma_struct = all_PFC_Contrast_gamma_struct_drug - all_PFC_Contrast_gamma_struct
left_PFC_Qmoney_gamma_struct = left_PFC_Contrast_gamma_struct_drug - left_PFC_Contrast_gamma_struct
right_PFC_Qmoney_gamma_struct = right_PFC_Contrast_gamma_struct_drug - right_PFC_Contrast_gamma_struct
mid_PFC_Qmoney_gamma_struct = mid_PFC_Contrast_gamma_struct_drug - mid_PFC_Contrast_gamma_struct