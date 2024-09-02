import os

def update_copyright(file_path, old_texts, new_text):
    with open(file_path,'r') as file:
        content = file.read()

    updated_content = content
    for old_text in old_texts:
        updated_content = updated_content.replace(old_text, new_text)

    if content != updated_content:
        with open(file_path, 'w') as file:
            file.write(updated_content)
        return True
    return False

def updated_files_From_list(file_list_path, old_texts, new_text):
    with open(file_list_path,'r') as file_list:
        files = file_list.readlines()

    changed_files_count = 0

    for file_path in files:
        file_path = file_path.strip()
        if os.path.isfile(file_path):
            if update_copyright(file_path,old_texts,new_text):
                changed_files_count += 1
                print(f"Updated: {file_path}")
            else:
                print(f"No changes made: {file_path}")

        else:
            print(f"File not found: {changed_files_count}")

    print(f"\nTotal files changed: {changed_files_count}")

if __name__ =="__main__":
    old_texts = ["Copyright 2020, General Electric Company","Copyright 2021, General Electric Company",
                 "Copyright 2022, General Electric Company","Copyright 2023, General Electric Company",
                 "Copyright 2023, GE HealthCare"
                 ]
    new_text = "Copyright 2024, GE Precision HealthCare"
    file_list_path = "./updateCopyrightTextFiles.txt"
    updated_files_From_list(file_list_path, old_texts, new_text)