#!/usr/bin/env python
# coding: utf-8

import argparse

def gfa2ctgpath(gfa):
    with open(gfa, "r") as f:
        gfa = {}
        for each in f:
            l = each.strip().split()
            if len(l) !=9 or l[0] != 'A' :
                continue
            if l[1] not in gfa:
                gfa[l[1]] = [each]
            else:
                gfa[l[1]].append(each)
    utg = {}
    for key in gfa.keys():
        utg[key] = []
    utg_len = dict.fromkeys(gfa.keys(),0)
    utg_overlap = dict.fromkeys(gfa.keys(),0)
    for key in gfa.keys():
        pre = None
        for l in gfa[key]:
            if pre == None :
                pre = l.strip().split()
                pre[2] = int(pre[2])
                pre[6] = int(pre[6])
                if  pre[3] == '-' :
                    utg[key].append(pre[4] + ':B')
                elif  pre[3] == '+' :
                    utg[key].append(pre[4] + ':E')
            else:
                l = l.strip().split()
                l[2] = int(l[2])
                l[6] = int(l[6])
                if  l[3] == '-' :
                    utg[key].append(l[4] + ':B')
                    utg_len[key] += l[2] - pre[2]
                    utg_overlap[key] += abs(l[6] - l[2] + pre[2] )
                elif  l[3] == '+' :
                    utg[key].append(l[4] + ':E')
                    utg_len[key] += l[2] - pre[2]
                    utg_overlap[key] += abs(l[6] - l[2] + pre[2] )
                pre = l
    with open("ctgpath","w") as f:
        for key in utg.keys():
            if(len(utg[key])>1):
                f.write("\t".join((key,utg[key][0],utg[key][1],utg[key][-1],str(utg_len[key]),str(utg_overlap[key]),"~".join(utg[key]))))
                f.write("\n")
            else:
                f.write("\t".join((key,utg[key][0],utg[key][0],utg[key][-1],str(utg_len[key]),str(utg_overlap[key]),"~".join(utg[key]))))
                f.write("\n")


def reverse_end(node_id):
    node_id, end = node_id.split(":")
    new_end = "B" if end == "E" else "E"
    return node_id + ":" + new_end


def find_chr_path(agp):
    chr_path = {}
    ctg_dir = {}
    with open(agp,'r') as f:
        for each in f:
            each = each.strip().split()
            if each[4] == "U":
                continue
            ctg_dir[each[5]] = each[8]
            if int(each[7]) <= 500000 :
                continue
            if each[0] not in chr_path:
                chr_path[each[0]] =  [each[5]]
            else:
                chr_path[each[0]].append(each[5])
    ctg_path = {} 
    with open("ctgpath",'r') as f:
        for each in f:
            each = each.strip().split()
            if each[0] not in ctg_dir.keys():
                continue
            if ctg_dir[each[0]] == '+':

                ctg_path[each[0]] = each
            else :
                path = [ reverse_end(i) for i in each[6].split("~")]
                path = list(reversed(path))
                s = path[0]
                if len(path) == 1:
                    v = path[0]
                else:
                    v = path[1]
                t = path[-1]
                path = '~'.join(path)
                ctg_path[each[0]] = [each[0],s,v,t,each[4],each[5],path]
    with open("ctg_paths","w") as f:
        for each  in ctg_path:
            each = ctg_path[each]
            f.write("\t".join(each))
            f.write("\n")
    f = open("chr_paths","w")
    for each  in chr_path.keys():
        name = each
        each = chr_path[name]
        begin = ctg_path[each[0]][1]
        v = ctg_path[each[0]][2]
        end = ctg_path[each[-1]][3]
        length = int(ctg_path[each[0]][4])
        score = int(ctg_path[each[0]][5])
        pre =  ctg_path[each[0]][3]
        path = []
        path.append(ctg_path[each[0]][6])
        #print each[0]
        i = 1
        while i < len(each):
            length += int(ctg_path[each[i]][4])
            score += int(ctg_path[each[i]][5])
            if pre == ctg_path[each[i]][1]:
                path.append(ctg_path[each[i]][6])
            else:
                path.append("~".join((pre,"gap",ctg_path[each[i]][1])))
                print((name,pre,"gap",ctg_path[each[i]][1]))
                path.append(ctg_path[each[i]][6])
            pre = ctg_path[each[i]][3]
            i += 1
        path = "|".join(path)
        f.write("\t".join((name,begin,v,end,str(length),str(score),path)))
        f.write("\n")
    f.close()


if __name__ == '__main__':
    parser =  argparse.ArgumentParser(usage = """chr_paths.py  --agp AGP --gfa GFA""",description = "get chr_paths using agp file and gfa file")
    parser.add_argument('--agp',type=str,help="agp file form contig anchoring",required = True)
    parser.add_argument('--gfa',type=str,help="gfa from hifiasm",required = True)
    arg = vars(parser.parse_args())
    agp = arg["agp"]
    gfa = arg["gfa"]

    gfa2ctgpath(gfa)
    find_chr_path(agp)
