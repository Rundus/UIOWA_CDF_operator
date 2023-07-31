# # --- scratchPaper.py ---
# # --- Author: C. Feltman ---
# # DESCRIPTION: script for general tinkering

# def read_file_with_space_delimiter(file_path):
#     try:
#         with open(file_path, 'r') as file:
#             # Read the file and split lines into a list of lines
#             lines = file.readlines()
#
#         data_list = []
#         for line in lines:
#             elements = line.strip().split()
#             # Convert the first three elements to integers
#             sublist = [int(elements[0]), int(elements[1]), int(elements[2]), str(elements[3]), int(elements[4])]
#             data_list.append(sublist)
#
#         return data_list
#     except FileNotFoundError:
#         print(f"Error: File '{file_path}' not found.")
#         return []
#
# if __name__ == "__main__":
#     file_path = r"C:\Users\cfelt\OneDrive\Desktop\testData.txt"  # Replace this with the actual path of your .txt file
#     data_list = read_file_with_space_delimiter(file_path)
#
#     # word_high, word_low, sign, '{0:016b}'.format(word_high)[-7:] + '{0:016b}'.format(word_low), fullword
#
#     for i in range(len(data_list)):
#         if data_list[i][2] == 1:
#             highWord = '{0:016b}'.format(data_list[i][0])[-7:]
#             lowWord = '{0:016b}'.format(data_list[i][1])
#             twoComp = highWord.replace('1', '2').replace('0', '1').replace('2', '0')
#             val = twoComp + lowWord
#
#             print(data_list[i][0],data_list[i][1],highWord, lowWord, twoComp, val,int(val,2))
#         else:
#             highWord ='{0:016b}'.format(data_list[i][0])[-7:]
#             lowWord = '{0:016b}'.format(data_list[i][1])
#             val =  highWord + lowWord
#             print(data_list[i][0],data_list[i][1],highWord,lowWord,val, int(val,2))

import numpy as np
a = [1,2,3]
print(np.delete(a,[0,2]))
