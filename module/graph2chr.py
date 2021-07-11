#!/usr/bin/python
import argparse
import logging
import sys
import networkx as nx
import re

RCMAP = dict(list(zip("ACGTacgtNn-", "TGCAtgcaNn-")))



def log(msg):
    sys.stderr.write(msg)
    sys.stderr.write('\n')


def rc(seq):
    return "".join([RCMAP[c] for c in seq[::-1]])

def reverse_end(node_id):
    node_id, end = node_id.split(":")
    new_end = "B" if end == "E" else "E"
    return node_id + ":" + new_end



def build_sg(sg_edges_list_fn, paf_fn, ctg_paths_fn):
    overlap_set = set()
    edges = []
    ctg_edges = {}
    with open(ctg_paths_fn,"r")  as f :
        for l in f:
            l = l.strip().split()
            edge = re.split("~|\|",l[6])
            edge = [ each.split(":")[0] for each in edge ]
            edges.extend(list(zip(edge[:-1], edge[1:])))
    for each in edges :
        overlap_pair = list(each)
        overlap_pair.sort()
        overlap_pair = tuple(overlap_pair)
        if overlap_pair not in overlap_set:
            overlap_set.add(overlap_pair)

    sg = nx.DiGraph()
    with open(paf_fn) as f:
        for line in f:
            l = line.strip().split()
            if len(l) != 12:
                continue
            f_id, g_id, score, identity = (l[5], l[0], int(l[7]) - int(l[8]), 99.9)
            g_s = 0 if l[4] == '+' else 1
            f_s, f_b, f_e, f_l = (int(i) for i in [0, l[7], l[8], l[6]])
            g_b, g_e, g_l = (int(i) for i in [l[2], l[3], l[1]])

            overlap_pair = [f_id, g_id]
            overlap_pair.sort()
            overlap_pair = tuple(overlap_pair)
            if overlap_pair not in overlap_set:  # don't allow duplicated records
                continue

            if g_s == 1:  # revered alignment, swapping the begin and end coordinates
                g_b, g_e = g_e, g_b

            # build the string graph edges for each overlap
            if f_b > 0:
                if g_b < g_e:
                    """
                         f.B         f.E
                      f  ----------->
                      g         ------------->
                                g.B           g.E
                    """
                    if f_b == 0 or g_e - g_l == 0:
                        continue
                    sg.add_edge("%s:B" % g_id, "%s:B" % f_id, label = (f_id, f_b, 0),
                                length=abs(f_b - 0),
                                e_score=-score,
                                identity=identity)
                    ctg_edges[("%s:B" % g_id, "%s:B" % f_id)] = ("%s:B" % g_id, "%s:B" % f_id,f_id, f_b, 0, -score, identity, "G")
                    sg.add_edge("%s:E" % f_id, "%s:E" % g_id, label = (g_id, g_e, g_l),
                                length=abs(g_e - g_l),
                                e_score=-score,
                                identity=identity)
                    ctg_edges[("%s:E" % f_id, "%s:E" % g_id)] = ("%s:E" % f_id, "%s:E" % g_id, g_id, g_e, g_l, -score, identity, "G")
                else:
                    """
                         f.B         f.E
                      f  ----------->
                      g         <-------------
                                g.E           g.B
                    """
                    if f_b == 0 or g_e == 0:
                        continue
                    sg.add_edge("%s:E" % g_id, "%s:B" % f_id, label = (f_id, f_b, 0),
                                length=abs(f_b - 0),
                                e_score=-score,
                                identity=identity)
                    ctg_edges[("%s:E" % g_id, "%s:B" % f_id)] = ("%s:E" % g_id, "%s:B" % f_id, f_id, f_b, 0, -score, identity, "G")
                    sg.add_edge("%s:E" % f_id, "%s:B" % g_id, label = (g_id, g_e, 0),
                                length=abs(g_e - 0),
                                e_score=-score,
                                identity=identity)
                    ctg_edges[("%s:E" % f_id, "%s:B" % g_id)] = ("%s:E" % f_id, "%s:B" % g_id, g_id, g_e, 0, -score, identity, "G")
            else:
                if g_b < g_e:
                    """
                                        f.B         f.E
                      f                 ----------->
                      g         ------------->
                                g.B           g.E
                    """
                    if g_b == 0 or f_e - f_l == 0:
                        continue
                    sg.add_edge("%s:B" % f_id, "%s:B" % g_id, label = (g_id, g_b, 0),
                                length=abs(g_b - 0),
                                e_score=-score,
                                identity=identity)
                    ctg_edges[("%s:B" % f_id, "%s:B" % g_id)] = ("%s:B" % f_id, "%s:B" % g_id, g_id, g_b, 0, -score, identity, "G")
                    sg.add_edge("%s:E" % g_id, "%s:E" % f_id, label = (f_id, f_e, f_l),
                                length=abs(f_e - f_l),
                                e_score=-score,
                                identity=identity)
                    ctg_edges[("%s:E" % g_id, "%s:E" % f_id)] = ("%s:E" % g_id, "%s:E" % f_id, f_id, f_e, f_l, -score, identity, "G")
                else:
                    """
                                        f.B         f.E
                      f                 ----------->
                      g         <-------------
                                g.E           g.B
                    """
                    if g_b - g_l == 0 or f_e - f_l == 0:
                        continue
                    sg.add_edge("%s:B" % f_id, "%s:E" % g_id, label = (g_id, g_b, g_l),
                                length=abs(g_b - g_l),
                                e_score=-score,
                                identity=identity)
                    ctg_edges[("%s:B" % f_id, "%s:E" % g_id)] = ("%s:B" % f_id, "%s:E" % g_id, g_id, g_b, g_l, -score, identity, "G")
                    sg.add_edge("%s:B" % g_id, "%s:E" % f_id, label = (f_id, f_e, f_l),
                                length=abs(f_e - f_l),
                                e_score=-score,
                                identity=identity)
                    ctg_edges[("%s:B" % g_id, "%s:E" % f_id)] = ("%s:B" % g_id, "%s:E" % f_id, f_id, f_e, f_l, -score, identity, "G")

    with open(sg_edges_list_fn) as f:
        for l in f:
            l = l.strip().split()
            """001039799:E 000333411:E 000333411 17524 20167 17524 99.62"""
            v, w, rid, s, t, aln_score, idt, type_ = l
            if  type_ == "TR":
                continue

            overlap_pair = [v, w]
            overlap_pair.sort()
            overlap_pair = tuple(overlap_pair)
            if overlap_pair  in overlap_set:  # don't allow duplicated records
                continue

            s = int(s)
            t = int(t)
            aln_score = int(aln_score)
            idt = float(idt)
            length = abs(s - t)
            label = (rid, s, t)
            sg.add_edge(v, w, label=label, length=length, e_score=aln_score, identity=idt)


    return sg,ctg_edges


def shortestpath(sg,ctg_id,start,t,edges,a_ctg_group,rd):
    all_alt_path = []
    c_graph = nx.ego_graph(sg, start, radius=rd)
    shortest_path = nx.shortest_path(c_graph, start, t, "e_score")
    length = len(shortest_path)
    score = nx.shortest_path_length(c_graph, start, t, "e_score")
    all_alt_path.append((length, shortest_path))
    count = 3
    while count > 0:
        n0 = shortest_path[0]
        for n1 in shortest_path[1:]:
            c_graph.remove_edge(n0, n1)
            n0 = n1
        try:
            shortest_path = nx.shortest_path(c_graph, start, t, "e_score")
            score = nx.shortest_path_length(c_graph, start, t, "e_score")
            all_alt_path.append((length, shortest_path))

        except nx.exception.NetworkXNoPath:
            break
        count -= 1
    all_alt_path.sort()
    all_alt_path.reverse()
    a_ctg_group[(ctg_id, start, t)] = all_alt_path

    for score, atig_path in all_alt_path:
        atig_path_edges = list(zip(atig_path[:-1], atig_path[1:]))
        for vv, ww in atig_path_edges:
            e = sg[vv][ww]
            rid, ss, tt = e["label"]
            aln_score = e["e_score"]
            idt = e["identity"]
            edges[(vv,ww)] = (vv, ww, rid, ss, tt, aln_score, idt, "G")

    shortest_path = all_alt_path[0][1]

    return shortest_path


def traverse_sg(ctg, sg):
    new_ctg = {}
    a_ctg_group = {}
    edges = {}
    for key in ctg.keys():
        l = ctg[key]
        ctg_id, s0, v0, t0, ctg_length, ctg_score, utgs = l
        ctg_length = int(ctg_length)
        ctg_score = int(ctg_score)
        utgs = utgs.split("|")
        utg_path = []
        tmp_utg_path = []

        l = utgs[0].split("~")
        utg_path.extend(l)
        i = 1
        while i < len(utgs):
            l = utgs[i].split("~")
            s = l[0]
            v = l[1]
            t0 = l[-1]
            t = l[-1]
            if v != "gap":
                utg_path.extend(l[1:])
            else:
                n_ego_graph = nx.ego_graph(sg, s, radius=100)
                s_path = []
                if t  in n_ego_graph.nodes():
                    try:
                        s_path = nx.shortest_path(n_ego_graph, s, t, "e_score")
                    except nx.exception.NetworkXNoPath:
                        s_path = []
                rd = 100
                if len(s_path) != 0:
                    gap_path = shortestpath(n_ego_graph,ctg_id, s, t, edges, a_ctg_group,rd)
                    utg_path.extend(gap_path[1:])
                else: 
                    tmp_utg_path = utg_path[:]

                    j = 1
                    while len(s_path) == 0 and i  + 1 < len(utgs) and j <=100  and len(tmp_utg_path) > 1:
                        m = utgs[i + 1].split("~")
                        if j >= len(m):
                            break
                        t = m[j]
                        tmp_utg_path.pop()
                        start = tmp_utg_path[-1]
                        if start == 'N':
                            break
                        n_ego_graph = nx.ego_graph(sg, start, radius=rd)
                        if t  in n_ego_graph.nodes():
                            try:
                                s_path = nx.shortest_path(n_ego_graph, start, t, "e_score")
                            except nx.exception.NetworkXNoPath:
                                s_path = []
                        if len(s_path) != 0:
                            print(ctg_id, start, t)
                            gap_path = shortestpath(n_ego_graph,ctg_id, start, t, edges, a_ctg_group, rd)
                            print(ctg_id, start, t,gap_path)
                            tmp_utg_path.extend(gap_path[1:])
                            tmp_utg_path.extend(m[j+1:])
                            i = i + 1
                            break
                        else:
                            j = j + 1
                        rd = rd + 7
                    if len(s_path) != 0:
                        utg_path = tmp_utg_path[:]
                    else :
                        utg_path.append('N')
                        print("N",ctg_id,s,v,t0)
                    

            i = i + 1

        new_ctg[ctg_id] = ctg_id, s0, v0, t, ctg_length, ctg_score, utg_path

    return new_ctg,a_ctg_group,edges



def yield_first_seq(one_path_edges, seqs):
    if one_path_edges and one_path_edges[0][0] != one_path_edges[-1][1]:
        (vv, ww) = one_path_edges[0]
        (vv_rid, vv_letter) = vv.split(":")
        if vv_letter == 'E':
            first_seq = seqs[vv_rid]
        else:
            assert vv_letter == 'B'
            first_seq = "".join([RCMAP[c] for c in seqs[vv_rid][::-1]])
        yield first_seq



def compose_ctg(seqs, edge_data, ctg_id, path_edges, proper_ctg):
    total_score = 0
    total_length = 0
    edge_lines = []
    sub_seqs = []

    if proper_ctg:
        sub_seqs = list(yield_first_seq(path_edges, seqs))
        total_length = 0 if len(sub_seqs) == 0 else len(sub_seqs[0])

    for vv, ww in path_edges:
        if vv == "N" or ww == "N":
            sub_seqs.append('N' * 50)
            edge_lines.append('%s %s %s %s %d %d %d %0.2f' % (
                ctg_id, vv, ww, vv, 0, 0, 0, 0))
            total_length += 50
        else:
            rid, s, t, aln_score, idt, e_seq = edge_data[(vv, ww)]
            sub_seqs.append(e_seq)
            edge_lines.append('%s %s %s %s %d %d %d %0.2f' % (
                ctg_id, vv, ww, rid, s, t, aln_score, idt))
            total_length += abs(s - t)
            total_score += aln_score

    return edge_lines, sub_seqs, total_score, total_length


def run(reads_fasta_fn, sg_edges_list_fn,paf_fn, ctg_paths_fn):
    sg,edges = build_sg(sg_edges_list_fn, paf_fn, ctg_paths_fn)


    ctg = {}
    with open(ctg_paths_fn) as f:
        for l in f:
            l = l.strip().split()
            ctg_id, s, v, t, ctg_length, ctg_score, utg = l
            ctg_id = ctg_id
            ctg[ctg_id] = ctg_id, s, v, t, ctg_length, ctg_score, utg
    ctg,gap_group,raw_edges = traverse_sg(ctg, sg)
    edges.update(raw_edges)
    with open("chr_path_new","w") as f:
        for key in ctg.keys():
            path = '~'.join(ctg[key][6])
            f.write('%s %s %s %s %d %d %s\n' % (ctg[key][0],ctg[key][1],ctg[key][2],ctg[key][3],ctg[key][4],ctg[key][5],path))

    reads_in_layout = set()
    for (v, w) in edges.keys():
        r1 = v.split(":")[0]
        reads_in_layout.add(r1)
        r2 = w.split(":")[0]
        reads_in_layout.add(r2)

    seqs = {}
    with open(reads_fasta_fn, 'r')  as f:
        for each in f:
            each = each.strip()
            if each[0] == '>':
                rn = each[1:]
            elif rn in reads_in_layout:
                seqs[rn] = each

    edge_data = {}
    for v, w in edges.keys():
        vv, ww, rid, ss, tt, aln_score, idt, type = edges[(v,w)]
        if ss < tt:
            e_seq = seqs[rid][ss:tt]
        else:
            e_seq = "".join([RCMAP[c] for c in seqs[rid][tt:ss][::-1]])
        edge_data[(v, w)] = (rid, ss, tt, aln_score, idt, e_seq)


    p_ctg_out = open("p_ctg.fasta", "w")
    p_ctg_t_out = open("p_ctg_tiling_path", "w")
    gap_out = open("gap_all.fasta", "w")
    gap_t_out = open("gap_all_tiling_path", "w")

    a_id = 0
    for (ctg_id,v, w) in gap_group.keys():
        atig_output = []
        for sub_id in range(len(gap_group[(ctg_id,v, w)])):
            score, atig_path = gap_group[(ctg_id,v, w)][sub_id]
            atig_path_edges = list(zip(atig_path[:-1], atig_path[1:]))

            a_ctg_id = '%s-%03d-%02d' % (ctg_id, a_id + 1, sub_id)
            a_edge_lines, sub_seqs, a_total_score, a_total_length = compose_ctg(
                seqs, edge_data, a_ctg_id, atig_path_edges, False)

            seq = ''.join(sub_seqs)
            atig_output.append((v, w, atig_path, a_total_length, a_total_score, seq, atig_path_edges, a_ctg_id,
                                a_edge_lines))

        for  data in atig_output:
            v, w, tig_path, a_total_length, a_total_score, seq, atig_path_edges, a_ctg_id, a_edge_lines = data

            gap_t_out.write('\n'.join(a_edge_lines))
            gap_t_out.write('\n')
            gap_out.write('>%s %s %s %d %d %d \n' % (
            a_ctg_id, v, w, a_total_length, a_total_score, len(atig_path_edges)))
            gap_out.write(''.join(seq))
            gap_out.write('\n')
        a_id += 1


    layout_ctg = set()
    for l in ctg.keys():
        ctg_id, s, v, t, ctg_length, ctg_score, utg = ctg[l]
        if (reverse_end(t), reverse_end(s)) in layout_ctg:
            continue
        else:
            layout_ctg.add((s, t))
        ctg_label = s + "~" + v + "~" + t
        length = int(ctg_length)

        total_score = ctg_score

        one_path_edges = list(zip(utg[:-1], utg[1:]))

        p_edge_lines, p_ctg_seq_chunks, p_total_score, p_total_length = compose_ctg(seqs, edge_data, ctg_id, one_path_edges, True)

        p_ctg_t_out.write('\n'.join(p_edge_lines))
        p_ctg_t_out.write('\n')

        p_ctg_out.write('>%s %s %d %d\n' % (ctg_id, ctg_label, p_total_length, total_score))
        p_ctg_out.write(''.join(p_ctg_seq_chunks))
        p_ctg_out.write('\n')

    p_ctg_out.close()
    p_ctg_t_out.close()
    gap_out.close()
    gap_t_out.close()


class HelpF(argparse.RawDescriptionHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    pass

def main(argv=sys.argv):
    description = 'Generate the primary and tiling paths, given the string graph.'
    epilog = """
We write these:

    p_ctg_out = open("p_ctg.fasta", "w")
    gap_out = open("gap_all.fasta", "w")
    p_ctg_t_out = open("p_ctg_tiling_path", "w")
    gap_t_out = open("gap_all_tiling_path", "w")
"""
    parser = argparse.ArgumentParser(
            description=description,
            formatter_class=HelpF,
            epilog=epilog)
    parser.add_argument('--reads-fasta-fn', type=str,
            help='Input. reads file generated by hifiasm.',required = True)
    parser.add_argument('--paf-fn', type=str,
                        help='Input. reads alignment file generated by hifiasm.',required = True)
    parser.add_argument('--sg-edges-list-fn', type=str,
            help='Input. File containing string graph edges, produced by ovlp2graph.py.',required = True)
    parser.add_argument('--ctg-paths-fn', type=str,
            help='Input. File containing contig paths.',required = True)
    args = parser.parse_args(argv[1:])
    print(vars(args))
    run(**vars(args))

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    main(sys.argv)
