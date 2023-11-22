import os


folder_data="../data/"

folder_norm_ds=folder_data+"NormDatasets/"


def erase_all_csvs(path):
    for file_path in os.listdir(path):
        if file_path.startswith("predictionsHybrid_"):
            os.remove(path+file_path)


erase_all_csvs(folder_norm_ds+"20000Prot/")
