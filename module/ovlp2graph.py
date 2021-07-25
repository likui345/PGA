import networkx as nx
import argparse
import logging
import os
import re
import sys
import collections.abc

# http://code.activestate.com/recipes/576694/
class OrderedSet(collections.abc.MutableSet):

    def __init__(self, iterable=None):
        self.end = end = []
        end += [None, end, end]         # sentinel node for doubly linked list
        self.map = {}                   # key --> [key, prev, next]
        if iterable is not None:
            self |= iterable

    def __len__(self):
        return len(self.map)

    def __contains__(self, key):
        return key in self.map

    def update(self, other):
        for i in other:
            self.add(i)

    def add(self, key):
        if key not in self.map:
            end = self.end
            curr = end[1]
            curr[2] = end[1] = self.map[key] = [key, curr, end]

    def discard(self, key):
        if key in self.map:
            key, prev, next = self.map.pop(key)
            prev[2] = next
            next[1] = prev

    def __iter__(self):
        end = self.end
        curr = end[2]
        while curr is not end:
            yield curr[0]
            curr = curr[2]

    def __reversed__(self):
        end = self.end
        curr = end[1]
        while curr is not end:
            yield curr[0]
            curr = curr[1]

    def pop(self, last=True):
        if not self:
            raise KeyError('set is empty')
        key = self.end[1][0] if last else self.end[2][0]
        self.discard(key)
        return key

    def __repr__(self):
        if not self:
            return '%s()' % (self.__class__.__name__,)
        return '%s(%r)' % (self.__class__.__name__, list(self))

    def __eq__(self, other):
        if isinstance(other, OrderedSet):
            return len(self) == len(other) and list(self) == list(other)
        return set(self) == set(other)



from collections import OrderedDict as dict
LOG = logging.getLogger(__name__)


class SGNode(object):
    """
    class representing a node in the string graph
    """

    def __init__(self, node_name):
        self.name = node_name
        self.out_edges = []
        self.in_edges = []

    def add_out_edge(self, out_edge):
        self.out_edges.append(out_edge)

    def add_in_edge(self, in_edge):
        self.in_edges.append(in_edge)


class SGEdge(object):
    """
    class representing an edge in the string graph
    """

    def __init__(self, in_node, out_node):
        self.in_node = in_node
        self.out_node = out_node
        self.attr = {}

    def set_attribute(self, attr, value):
        self.attr[attr] = value


def reverse_end(node_name):
    if (node_name == 'NA'):
        return node_name
    if (len(node_name) < 2 or (node_name[-2:] not in [':B', ':E'])):
        raise Exception(
            'Invalid node name. Node name passed to method: "{node_name}", expected format: "(%d)+:[BE]" or "NA".'.format(node_name=node_name))
    node_id, end = node_name.split(":")
    new_end = "B" if end == "E" else "E"
    return node_id + ":" + new_end


class StringGraph(object):
    """
    class representing the string graph
    """

    def __init__(self):
        self.nodes = {}
        self.edges = {}
        self.e_reduce = {}
        self.best_in = {}

    def add_node(self, node_name):
        """
        add a node into the graph by given a node name
        """
        if node_name not in self.nodes:
            self.nodes[node_name] = SGNode(node_name)

    def add_edge(self, in_node_name, out_node_name, **attributes):
        """
        add an edge into the graph by given a pair of nodes
        """
        if (in_node_name, out_node_name) not in self.edges:

            self.add_node(in_node_name)
            self.add_node(out_node_name)
            in_node = self.nodes[in_node_name]
            out_node = self.nodes[out_node_name]

            edge = SGEdge(in_node, out_node)
            self.edges[(in_node_name, out_node_name)] = edge
            in_node.add_out_edge(edge)
            out_node.add_in_edge(edge)
        edge = self.edges[(in_node_name, out_node_name)]
        for (k, v) in attributes.items():
            edge.attr[k] = v

    def init_reduce_dict(self):
        for e in self.edges:
            self.e_reduce[e] = False

    def bfs_nodes(self, n, exclude=None, depth=5):
        all_nodes = set()
        all_nodes.add(n)
        candidate_nodes = set()
        candidate_nodes.add(n)
        dp = 1
        while dp < depth and len(candidate_nodes) > 0:
            v = candidate_nodes.pop()
            for e in v.out_edges:
                w = e.out_node
                if w == exclude:
                    continue
                if w not in all_nodes:
                    all_nodes.add(w)
                    if len(w.out_edges) > 0:
                        candidate_nodes.add(w)
            dp += 1

        return all_nodes

    def mark_chimer_edges(self):

        multi_in_nodes = {}
        multi_out_nodes = {}
        for n_name in self.nodes:
            n = self.nodes[n_name]
            out_nodes = [e.out_node for e in n.out_edges if self.e_reduce[(
                e.in_node.name, e.out_node.name)] == False]
            in_nodes = [e.in_node for e in n.in_edges if self.e_reduce[(
                e.in_node.name, e.out_node.name)] == False]

            if len(out_nodes) >= 2:
                multi_out_nodes[n_name] = out_nodes
            if len(in_nodes) >= 2:
                multi_in_nodes[n_name] = in_nodes

        chimer_candidates = set()
        out_set = set()
        in_set = set()
        for n_name in multi_out_nodes:
            out_nodes = set(multi_out_nodes[n_name])
            out_set |= out_nodes

        for n_name in multi_in_nodes:
            in_nodes = set(multi_in_nodes[n_name])
            in_set |= in_nodes

        chimer_candidates = out_set & in_set

        chimer_nodes = []
        chimer_edges = set()
        for n in chimer_candidates: # sort, or OrderedSet
            out_nodes = set([e.out_node for e in n.out_edges])
            test_set = set()
            for in_node in [e.in_node for e in n.in_edges]:
                test_set = test_set | set(
                    [e.out_node for e in in_node.out_edges])
            test_set -= set([n])
            if len(out_nodes & test_set) == 0:
                flow_node1 = set()
                flow_node2 = set()
                for v in list(out_nodes):
                    flow_node1 |= self.bfs_nodes(v, exclude=n)
                for v in list(test_set):
                    flow_node2 |= self.bfs_nodes(v, exclude=n)
                if len(flow_node1 & flow_node2) == 0:
                    for e in n.out_edges:
                        v, w = e.in_node.name, e.out_node.name
                        if self.e_reduce[(v, w)] != True:
                            self.e_reduce[(v, w)] = True
                            chimer_edges.add((v, w))
                            rv = reverse_end(w)
                            rw = reverse_end(v)
                            self.e_reduce[(rv, rw)] = True
                            chimer_edges.add((rv, rw))

                    for e in n.in_edges:
                        v, w = e.in_node.name, e.out_node.name
                        if self.e_reduce[(v, w)] != True:
                            self.e_reduce[(v, w)] = True
                            chimer_edges.add((v, w))
                            rv = reverse_end(w)
                            rw = reverse_end(v)
                            self.e_reduce[(rv, rw)] = True
                            chimer_edges.add((rv, rw))
                    chimer_nodes.append(n.name)
                    chimer_nodes.append(reverse_end(n.name))

        return chimer_nodes, chimer_edges

    def mark_spur_edge(self):

        removed_edges = set()
        for v in self.nodes:
            if len([e for e in self.nodes[v].out_edges if self.e_reduce[(e.in_node.name, e.out_node.name)] != True]) > 1:
                for out_edge in self.nodes[v].out_edges:
                    w = out_edge.out_node.name

                    if len(self.nodes[w].out_edges) == 0 and self.e_reduce[(v, w)] != True:
                        self.e_reduce[(v, w)] = True
                        removed_edges.add((v, w))
                        v2, w2 = reverse_end(w), reverse_end(v)
                        self.e_reduce[(v2, w2)] = True
                        removed_edges.add((v2, w2))

            if len([e for e in self.nodes[v].in_edges if self.e_reduce[(e.in_node.name, e.out_node.name)] != True]) > 1:
                for in_edge in self.nodes[v].in_edges:
                    w = in_edge.in_node.name
                    if len(self.nodes[w].in_edges) == 0 and self.e_reduce[(w, v)] != True:
                        self.e_reduce[(w, v)] = True
                        removed_edges.add((w, v))
                        v2, w2 = reverse_end(w), reverse_end(v)
                        self.e_reduce[(w2, v2)] = True
                        removed_edges.add((w2, v2))
        return removed_edges

    def mark_tr_edges(self):
        """
        transitive reduction
        """
        n_mark = {}
        e_reduce = self.e_reduce
        FUZZ = 500
        for n in self.nodes:
            n_mark[n] = "vacant"

        for (n_name, node) in self.nodes.items():

            out_edges = node.out_edges
            if len(out_edges) == 0:
                continue

            out_edges.sort(key=lambda x: x.attr["length"])

            for e in out_edges:
                w = e.out_node
                n_mark[w.name] = "inplay"

            max_len = out_edges[-1].attr["length"]

            max_len += FUZZ

            for e in out_edges:
                e_len = e.attr["length"]
                w = e.out_node
                if n_mark[w.name] == "inplay":
                    w.out_edges.sort(key=lambda x: x.attr["length"])
                    for e2 in w.out_edges:
                        if e2.attr["length"] + e_len < max_len:
                            x = e2.out_node
                            if n_mark[x.name] == "inplay":
                                n_mark[x.name] = "eliminated"

            for e in out_edges:
                e_len = e.attr["length"]
                w = e.out_node
                w.out_edges.sort(key=lambda x: x.attr["length"])
                if len(w.out_edges) > 0:
                    x = w.out_edges[0].out_node
                    if n_mark[x.name] == "inplay":
                        n_mark[x.name] = "eliminated"
                for e2 in w.out_edges:
                    if e2.attr["length"] < FUZZ:
                        x = e2.out_node
                        if n_mark[x.name] == "inplay":
                            n_mark[x.name] = "eliminated"

            for out_edge in out_edges:
                v = out_edge.in_node
                w = out_edge.out_node
                if n_mark[w.name] == "eliminated":
                    e_reduce[(v.name, w.name)] = True
                    v_name, w_name = reverse_end(w.name), reverse_end(v.name)
                    e_reduce[(v_name, w_name)] = True
                n_mark[w.name] = "vacant"

    def mark_best_overlap(self):
        """
        find the best overlapped edges
        """

        best_edges = set()
        removed_edges = set()

        for v in self.nodes:

            out_edges = self.nodes[v].out_edges
            if len(out_edges) > 0:
                out_edges.sort(key=lambda e: -e.attr["score"])
                for e in out_edges:
                    if self.e_reduce[(e.in_node.name, e.out_node.name)] != True:
                        best_edges.add((e.in_node.name, e.out_node.name))
                        break

            in_edges = self.nodes[v].in_edges
            if len(in_edges) > 0:
                in_edges.sort(key=lambda e: -e.attr["score"])
                for e in in_edges:
                    if self.e_reduce[(e.in_node.name, e.out_node.name)] != True:
                        best_edges.add((e.in_node.name, e.out_node.name))
                        self.best_in[v] = e.in_node.name
                        break

        LOG.debug(f"X {len(best_edges)}")

        for (e_n, e) in self.edges.items():
            v = e_n[0]
            w = e_n[1]
            if self.e_reduce[(v, w)] != True:
                if (v, w) not in best_edges:
                    self.e_reduce[(v, w)] = True
                    removed_edges.add((v, w))
                    v2, w2 = reverse_end(w), reverse_end(v)
                    self.e_reduce[(v2, w2)] = True
                    removed_edges.add((v2, w2))

        return removed_edges

    def resolve_repeat_edges(self):

        edges_to_reduce = []
        nodes_to_test = set()
        for (v_n, v) in self.nodes.items():

            out_nodes = []
            for e in v.out_edges:
                if self.e_reduce[(e.in_node.name, e.out_node.name)] == False:
                    out_nodes.append(e.out_node.name)

            in_nodes = []
            for e in v.in_edges:
                if self.e_reduce[(e.in_node.name, e.out_node.name)] == False:
                    in_nodes.append(e.in_node.name)

            if len(out_nodes) == 1 and len(in_nodes) == 1:
                nodes_to_test.add(v_n)

        for v_n in list(nodes_to_test):

            v = self.nodes[v_n]

            out_nodes = []
            for e in v.out_edges:
                if self.e_reduce[(e.in_node.name, e.out_node.name)] == False:
                    out_nodes.append(e.out_node.name)

            in_nodes = []
            for e in v.in_edges:
                if self.e_reduce[(e.in_node.name, e.out_node.name)] == False:
                    in_nodes.append(e.in_node.name)

            in_node_name = in_nodes[0]

            for out_edge in self.nodes[in_node_name].out_edges:
                vv = out_edge.in_node.name
                ww = out_edge.out_node.name

                ww_out = self.nodes[ww].out_edges
                v_out = self.nodes[v_n].out_edges
                ww_out_nodes = set([n.out_node.name for n in ww_out])
                v_out_nodes = set([n.out_node.name for n in v_out])
                o_overlap = len(ww_out_nodes & v_out_nodes)

                ww_in_count = 0
                for e in self.nodes[ww].in_edges:
                    if self.e_reduce[(e.in_node.name, e.out_node.name)] == False:
                        ww_in_count += 1

                if ww != v_n and\
                   self.e_reduce[(vv, ww)] == False and\
                   ww_in_count > 1 and\
                   ww not in nodes_to_test and\
                   o_overlap == 0:
                    edges_to_reduce.append((vv, ww))

            out_node_name = out_nodes[0]

            for in_edge in self.nodes[out_node_name].in_edges:
                vv = in_edge.in_node.name
                ww = in_edge.out_node.name

                vv_in = self.nodes[vv].in_edges
                v_in = self.nodes[v_n].in_edges
                vv_in_nodes = set([n.in_node.name for n in vv_in])
                v_in_nodes = set([n.in_node.name for n in v_in])
                i_overlap = len(vv_in_nodes & v_in_nodes)

                vv_out_count = 0
                for e in self.nodes[vv].out_edges:
                    if self.e_reduce[(e.in_node.name, e.out_node.name)] == False:
                        vv_out_count += 1

                if vv != v_n and\
                   self.e_reduce[(vv, ww)] == False and\
                   vv_out_count > 1 and\
                   vv not in nodes_to_test and\
                   i_overlap == 0:
                    edges_to_reduce.append((vv, ww))

        removed_edges = set()
        for e in edges_to_reduce:
            self.e_reduce[e] = True
            removed_edges.add(e)

        return removed_edges

def reverse_edge(e):
    e1, e2 = e
    return reverse_end(e2), reverse_end(e1)


def reverse_path(p):
    p = p[::-1]
    return [reverse_end(n) for n in p]



def init_string_graph(overlap_data):
    sg = StringGraph()

    overlap_set = set()
    for od in overlap_data:
        f_id, g_id, score, identity = od[:4]
        f_s, f_b, f_e, f_l = od[4:8]
        g_s, g_b, g_e, g_l = od[8:12]
        overlap_pair = [f_id, g_id]
        overlap_pair.sort()
        overlap_pair = tuple(overlap_pair)
        if overlap_pair in overlap_set:  # don't allow duplicated records
            continue
        else:
            overlap_set.add(overlap_pair)

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
                sg.add_edge("%s:B" % g_id, "%s:B" % f_id, label="%s:%d-%d"%(f_id, f_b, 0),
                            length=abs(f_b - 0),
                            score=-score,
                            identity=identity)
                sg.add_edge("%s:E" % f_id, "%s:E" % g_id, label="%s:%d-%d"%(g_id, g_e, g_l),
                            length=abs(g_e - g_l),
                            score=-score,
                            identity=identity)
            else:
                """
                     f.B         f.E
                  f  ----------->
                  g         <-------------
                            g.E           g.B
                """
                if f_b == 0 or g_e == 0:
                    continue
                sg.add_edge("%s:E" % g_id, "%s:B" % f_id, label="%s:%d-%d"%(f_id, f_b, 0),
                            length=abs(f_b - 0),
                            score=-score,
                            identity=identity)
                sg.add_edge("%s:E" % f_id, "%s:B" % g_id, label="%s:%d-%d"%(g_id, g_e, 0),
                            length=abs(g_e - 0),
                            score=-score,
                            identity=identity)
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
                sg.add_edge("%s:B" % f_id, "%s:B" % g_id, label="%s:%d-%d"%(g_id, g_b, 0),
                            length=abs(g_b - 0),
                            score=-score,
                            identity=identity)
                sg.add_edge("%s:E" % g_id, "%s:E" % f_id, label="%s:%d-%d"%(f_id, f_e, f_l),
                            length=abs(f_e - f_l),
                            score=-score,
                            identity=identity)
            else:
                """
                                    f.B         f.E
                  f                 ----------->
                  g         <-------------
                            g.E           g.B
                """
                if g_b - g_l == 0 or f_e - f_l == 0:
                    continue
                sg.add_edge("%s:B" % f_id, "%s:E" % g_id, label="%s:%d-%d"%(g_id, g_b, g_l),
                            length=abs(g_b - g_l),
                            score=-score,
                            identity=identity)
                sg.add_edge("%s:B" % g_id, "%s:E" % f_id, label="%s:%d-%d"%(f_id, f_e, f_l),
                            length=abs(f_e - f_l),
                            score=-score,
                            identity=identity)

    sg.init_reduce_dict()
    sg.mark_tr_edges()  # mark those edges that transitive redundant
    return sg

re_label = re.compile(r"(.*):(\d+)-(\d+)")


def init_digraph(sg, chimer_edges, removed_edges, spur_edges):
    with open("sg_edges_list", "w") as out_f:
        for v, w in sg.edges: # sort, or OrderedDict
            e = sg.edges[(v, w)]
            label = e.attr["label"]
            score = e.attr["score"]
            identity = e.attr["identity"]
            length = e.attr["length"]
            try:
                mo = re_label.search(label)
                rid = mo.group(1)
                sp = int(mo.group(2))
                tp = int(mo.group(3))
            except Exception:
                msg = 'parsing label="{}"'.format(label)
                LOG.exception(msg)
                raise
            assert length == abs(sp - tp)

            if not sg.e_reduce[(v, w)]:
                type_ = "G"
            elif (v, w) in chimer_edges:
                type_ = "C"
            elif (v, w) in removed_edges:
                type_ = "R"
            elif (v, w) in spur_edges:
                type_ = "S"
            else:
                assert sg.e_reduce[(v, w)]
                type_ = "TR"

            line = '%s %s %s %5d %5d %5d %5.2f %s' % (
                v, w, rid, sp, tp, score, identity, type_)
            print(line, file=out_f)



def yield_from_overlap_file(overlap_file):

    with open(overlap_file) as f:
        for line in f:
            if line.startswith('-'):
                break
            l = line.strip().split()
            if len(l) != 12 :
                continue
            f_id, g_id, score, identity =  ( l[5],l[0],int(l[7]) - int(l[8]),99.9 )
            if score >= -500:
                continue
            g_strand = 0 if l[4] == '+' else 1
            f_strand, f_start, f_end, f_len = ( int(i) for i  in [ 0, l[7], l[8], l[6]])
            g_start, g_end, g_len = ( int(i) for i in [l[2], l[3], l[1]] )
            yield (f_id, g_id, score, identity,
                                f_strand, f_start, f_end, f_len,
                                g_strand, g_start, g_end, g_len)

def generate_nx_string_graph(sg, lfc=False, disable_chimer_bridge_removal=False):
    LOG.debug("{}".format(sum([1 for c in sg.e_reduce.values() if c])))
    LOG.debug("{}".format(sum([1 for c in sg.e_reduce.values() if not c])))

    if not disable_chimer_bridge_removal:
        chimer_nodes, chimer_edges = sg.mark_chimer_edges()

        with open("chimers_nodes", "w") as f:
            for n in chimer_nodes:
                print(n, file=f)
        del chimer_nodes
    else:
        chimer_edges = set()  # empty set

    spur_edges = sg.mark_spur_edge()

    removed_edges = set()
    if lfc == True:
        removed_edges = sg.resolve_repeat_edges()
    else:
        # mark those edges that are best overlap edges
        removed_edges = sg.mark_best_overlap()

    spur_edges.update(sg.mark_spur_edge())

    LOG.debug('{}'.format(sum([1 for c in sg.e_reduce.values() if not c])))

    init_digraph(sg, chimer_edges, removed_edges, spur_edges)




def ovlp_to_graph(args):
    overlap_data = yield_from_overlap_file(args.overlap_file)
    sg = init_string_graph(overlap_data)

    # remove spurs, remove putative edges caused by repeats
    generate_nx_string_graph(sg)


class HelpF(argparse.RawTextHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    pass

def main(argv=sys.argv):
    epilog = """
Outputs:
    - sg_edges_list
    - chimer_nodes 
"""
    parser = argparse.ArgumentParser(
            description='string graph assembler that is desinged for handling diploid genomes',
            epilog=epilog,
            formatter_class=HelpF)
    parser.add_argument(
        '--overlap-file', help='the filtered overlap data from step2.',required = True)

    args = parser.parse_args(argv[1:])
    logging.basicConfig(level=logging.INFO, stream=sys.stdout, format='%(msg)s')
    ovlp_to_graph(args)


if __name__ == "__main__":
    main(sys.argv)
