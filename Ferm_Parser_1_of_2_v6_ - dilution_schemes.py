import argparse
import csv
import math
import string
import os
from os import listdir
#from os import listdir, isfile, join, path

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-sp", type=str, required=True)
    parser.add_argument("-std", type=str, required=True)
    parser.add_argument("-d", type=str, required=True)
    parser.add_argument("-ds", type=str, required=True)
    parser.add_argument("-psd", type=str, required=True)
    parser.add_argument("-i", type=str, required=True)
    args = parser.parse_args()
    return args

def get_sp_dict(sp_dir):
    sp_dict = {}
    sp_paths = [os.path.join(sp_dir, f) for f in listdir(sp_dir)]
    for f in sp_paths:
        SP_ID = os.path.basename(f)[:-4]
        sp_dict[SP_ID] = {}
        reader = csv.reader(open(f, "r"), delimiter="\t")
        header = next(reader)
        for line in reader:
            row = line[0]
            sp_dict[SP_ID][row] = {}
            col = 1
            while col < len(line):
                ID = line[col]
                sp_dict[SP_ID][row][col] = ID
                col += 1
    return sp_dict

def get_standards_dict(std_file):
    std_dict = {}
    reader = csv.reader(open(std_file, "r"), delimiter="\t")
    header = next(reader)
    for line in reader:
        key = line[0]
        value = line[1]
        std_dict[key] = value
    return std_dict

def get_dilutions_dict(dil_dir):
    dil_dict = {}
    dil_paths = [os.path.join(dil_dir, f) for f in listdir(dil_dir)]
    for f in dil_paths:
        dil_ID = os.path.basename(f)[:-4].split("_")[-1]
        dil_dict[dil_ID] = {}
        reader = csv.reader(open(f, "r"), delimiter="\t")
        header = next(reader)
        for line in reader:
            dil_count = line[0]
            dil_dict[dil_ID][dil_count] = line[1:]
    return dil_dict

def get_dilution_scheme_dict(ds_dir):
    """
    The dilution scheme dictionary will be separated into two text files. Therefore, this function with adjusted to
    pull in a directory containing both files. The dictionary built from the data of each file will contain
    a scheme as the main key. 

    """
    ds_dict = {}
    ds_paths = [os.path.join(ds_dir, f) for f in listdir(ds_dir)]
    ds_scheme = 1
    for ds_file in ds_paths:
        ds_dict[ds_scheme] = {}
        reader = csv.reader(open(ds_file, "r"), delimiter="\t")
        header = next(reader)
        for line in reader:
            row = line[0]
            ds_dict[ds_scheme][row] = {}
            col = 1
            while col < len(line):
                info_list = line[col].split(".")
                sr = info_list[0][0]
                sc = int(info_list[0][1:])
                dil = info_list[1][0]

                info = [sr, sc, dil]

                ds_dict[ds_scheme][row][col] = info
                col += 1
        ds_scheme += 1
    return ds_dict

def get_plate_scanner_dict(psd_file):
    ps_dict = {}
    reader = csv.reader(open(psd_file, "r") , delimiter="\t")
    header = next(reader)
    for line in reader:
        sp = line[0]
        dil = line[1]
        plates = line[2].split(",")
        for plate in plates:
            ps_dict[plate] = [sp, dil]
    return ps_dict

def get_data_dict(file):
    
    reader = csv.reader(open(file, "r", encoding="cp850"), delimiter=",")
    data_dict = {}
    a_D_ticker = 0
    a_D_key = "NA"
    # print(len(reader))
    while True:
        try:
            line = next(reader)
        except StopIteration:
            break

        if len(line) < 1: continue
        if line[0] != "Plate" and line[0] != "A": continue
        if line[0] == "Plate" and line[1] != "": #Keep in mind
            line = next(reader)
            plate = "P"+line[0]
            if plate not in data_dict: data_dict[plate] = {}

            a_D_ticker += 1
            if a_D_ticker % 2 == 1:
                a_D_key = "alpha"
            else:
                a_D_key = "DNA"
            data_dict[plate][a_D_key] = {}
            continue
        if line[0] == "A" and line[1] != "- ": #Keep in mind
            while len(line) > 1:
                row = line[0]
                data_dict[plate][a_D_key][row] = {}
                col = 1
                while col < 25:
                    well_value = line[col]
                    data_dict[plate][a_D_key][row][col] = well_value
                    col += 1
                    
                line = next(reader)
    return data_dict

def main(file):
    """
    This section had been adjusted to accommodate a new dilution scheme completed in the lab. Instead of source plates half
    filled with samples, with 8 point dilution occurring for samples on one 384 well greiner plate, source plates will be
    filled entirely with samples (and standard), where 8 point dilution will occur across two separate 384 well greiner
    plates for each sample. The scheme to use will be determined by whether the plate is odd or even numbered. 

    """
    sp_dict = get_sp_dict(args.sp)
    std_dict = get_standards_dict(args.std)
    dil_dict = get_dilutions_dict(args.d)
    ds_dict = get_dilution_scheme_dict(args.ds)
    ps_dict = get_plate_scanner_dict(args.psd)

    data_dict = get_data_dict(file)

    outname = file[:-4] + "_Reformat.txt"
    outhandle = open(outname, "w")
    header = ["sample_ID", "ferm_run", "sample_type", "nM", "BRno", "time", "pellet", "plate",
              "splateID", "splate row", "splate col", "row", "col", "dilution#", "volume", "alpha", "DNA"]
    outhandle.write("\t".join(header) + "\n")

    for plate in data_dict:
        if int(plate[1:]) % 2 == 0: 
            ds_scheme = 2
        else:
            ds_scheme = 1
        if plate == "P5": break
        SP_ID = ps_dict[plate][0]
        if SP_ID == "dummy": continue
        dil_ID = ps_dict[plate][1]
        for row in data_dict[plate]["alpha"]:
            for col in data_dict[plate]["alpha"][row]:
                alpha = data_dict[plate]["alpha"][row][col]
                DNA = data_dict[plate]["DNA"][row][col]
                ds_info = ds_dict[ds_scheme][row][col]
                sr = ds_info[0]
                sc = ds_info[1]
                dil = ds_info[2]
                dil_volume = dil_dict[dil_ID][dil][1]

                sample_ID = sp_dict[SP_ID][sr][sc]
                if sample_ID == "0": continue
                if "Stnd" not in sample_ID:
                    sample_info = sample_ID.split("_")
                    ferm_run = sample_info[0]
                    sample_type = "Experimental"
                    nM = "NA"
                    BRno = sample_info[1]
                    time = sample_info[2]
                    pellet = sample_info[-1]
                if "Stnd" in sample_ID:
                    sample_info = sample_ID.split("_")
                    ferm_run = sample_info[0]
                    sample_type = "Standard"
                    nM = std_dict[sample_info[0]]
                    BRno = "NA"
                    time = "NA"
                    pellet = "NA"

                outlist = [sample_ID, ferm_run, sample_type, nM, BRno, time, pellet, plate,
                           sr+str(sc), sr, sc, row, col, dil, dil_volume, alpha, DNA]
                outline = "\t".join([str(s) for s in outlist])
                outhandle.write(outline + "\n")

def main_2(file):
    """
    In order to sort 8 point dilution by sample, dictionary was adjusted to set keys by sample, not plate. 
    """
    file_2 = file[:-4] + "_Reformat.txt"

    reader = csv.reader(open(file_2, "r"), delimiter="\t")
    data_dict = {}
    header = next(reader)
    for line in reader:
        plate = line[7]
        
        sample_ID = line[0]
        dilution = int(line[13])
        # if plate not in data_dict: data_dict[plate] = {}
        # if sample_ID not in data_dict[plate]: data_dict[plate][sample_ID] = []

        if sample_ID not in data_dict: data_dict[sample_ID] = []
        temp_line = [dilution] + line
        # data_dict[plate][sample_ID].append(temp_line)
        data_dict[sample_ID].append(temp_line)

    outname_2 = file_2[:-4] + "_2.txt"
    outhandle_2 = open(outname_2, "w")
    outhandle_2.write("\t".join(header) + "\n")
    # for plate in data_dict:
    # for sample in data_dict[plate]:
    for sample in data_dict:
        # data_dict[plate][sample].sort()
        data_dict[sample].sort()

        # for line in data_dict[plate][sample]:
        for line in data_dict[sample]:
            outlist = line[1:]
            outline = "\t".join(outlist)
            outhandle_2.write(outline + "\n")        

    outname_3 = file_2[:-4] + "_3.txt"
    outhandle_3 = open(outname_3, "w")
    outheader = header + ["alpha_slope", "DNA_slope", "max_alpha_slope", "max_DNA_slope"]
    outhandle_3.write("\t".join(outheader) + "\n")
    # for plate in data_dict:
    # for sample in data_dict[plate]:
    #     working_data = data_dict[plate][sample]
    
    # Define the range to select max slope
    no_range = True
    while no_range:
        try:
            selected_range = int(input("Select range: "))
            no_range = False
        except ValueError:
            print("Invalid range!")
            
    for sample in data_dict:
        working_data = data_dict[sample]
        len_sample = len(working_data)
        # print(working_data)
        # for line in working_data:
        #     print(line)
        # print(len_sample)
        # if len_sample == 8:
            # print(working_data)
        i = 0
        j = 1
        while j < len_sample:
            vol_1 = float(working_data[i][-3])
            vol_2 = float(working_data[j][-3])
            alpha_1 = float(working_data[i][-2])
            alpha_2 = float(working_data[j][-2])
            DNA_1 = float(working_data[i][-1])
            DNA_2 = float(working_data[j][-1])
            try:
                alpha_slope = (alpha_2-alpha_1)/(vol_2-vol_1)
                DNA_slope = (DNA_2-DNA_1)/(vol_2-vol_1)
            except ZeroDivisionError:
                pass

            # data_dict[plate][sample][i].extend([alpha_slope, DNA_slope])
            data_dict[sample][i].extend([alpha_slope, DNA_slope])
            # print(data_dict[plate][sample][i])
            # print(data_dict[plate][sample][i])
            i += 1
            j += 1
        # print(max([float(data_dict[plate][sample][x][-2]) for x in range(4)]))
        # max_alpha_slope = max([data_dict[plate][sample][x][-2] for x in range(4)]) #I was told to use 4, but maybe that should change
        # max_DNA_slope = max([data_dict[plate][sample][x][-1] for x in range(4)])
        # data_dict[plate][sample][0].extend([max_alpha_slope, max_DNA_slope])       

        max_alpha_slope = max([data_dict[sample][x][-2] for x in range(selected_range)]) #I was told to use 4, but maybe that should change
        max_DNA_slope = max([data_dict[sample][x][-1] for x in range(selected_range)])
        data_dict[sample][0].extend([max_alpha_slope, max_DNA_slope])

        # for line in data_dict[plate][sample]:
        for line in data_dict[sample]:
            outlist = [str(s) for s in line[1:]]
            outline = "\t".join(outlist)
            outhandle_3.write(outline + "\n")

    outname_exp = file_2[:-4] + "_exp.txt"
    outname_std = file_2[:-4] + "_std.txt"
    outhandle_exp = open(outname_exp, "w")
    outhandle_std = open(outname_std, "w")
    outhandle_exp.write("\t".join(outheader) + "\n")
    outhandle_std.write("\t".join(outheader) + "\n")
    # for plate in data_dict:
    # for sample in data_dict[plate]:
    for sample in data_dict:
        # for line in data_dict[plate][sample]:
        for line in data_dict[sample]:
            sample_type = line[3]
            outlist = [str(s) for s in line[1:]]
            outline = "\t".join(outlist)
            if sample_type == "Experimental":
                outhandle_exp.write(outline + "\n")        
            if sample_type == "Standard":
                outhandle_std.write(outline + "\n")
            break


if __name__ == "__main__":
    args = get_args()
    file = args.i
    main(file)
    main_2(file)
