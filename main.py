import os

def remove_txt_files(folder_path):
    for file in os.listdir(folder_path):
        if file.endswith(".txt"):
            os.remove(os.path.join(folder_path, file))

remove_txt_files(os.path.join(os.getcwd(), "Logs"))