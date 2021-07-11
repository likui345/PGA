#!/usr/bin/python
import argparse

def filter(overlap_file, bestn, output):
    overlap_data = []
    overlap = {}
    fo = open(output,"w")
    with open(overlap_file) as f:
        for line in f:
            l = line.strip().split()
            if len(l) != 12 :
                continue
            score = int(l[8]) - int(l[7])
            if score < 2000:
                continue
            if int(l[1]) < 8000 or int(l[6]) < 8000 or int(l[1]) == int(l[3]) - int(l[2]) or int(l[6]) == int(l[8]) - int(l[7]):
                continue
            if l[0] not in overlap:
                if len(overlap_data) != 0:
                    for each in overlap_data:
                        fo.write('\t'.join(each))
                        fo.write('\n')
                overlap_data = [tuple(l)]
                overlap[l[0]] = 1
            else:
                overlap_data.append(tuple(l))
                if len(overlap_data) > bestn:
                    overlap_data.sort(key=lambda x: int(x[2]) - int(x[3]))
                    overlap_data = overlap_data[0:bestn]




# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    parser =  argparse.ArgumentParser(usage = """
    filter.py  --paf_fn PAF --bestn NUM --output OUTPUT
    """,description = "overlap data filter")
    parser.add_argument('--paf_fn',type=str,help="Input. reads alignment file",required = True)
    parser.add_argument('--bestn',type=int,default=20,help="output at least best n overlaps on 5' or 3' ends if possible.")
    parser.add_argument('--output',type=str,help="Output filename",required = True)
    arg = vars(parser.parse_args())
    overlap_file = arg["paf_fn"]
    bestn = arg["bestn"]
    output = arg["output"]
    filter(overlap_file,bestn,output)
