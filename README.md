{
 "cells": [
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "794673f1",
   "metadata": {},
   "outputs": [],
   "source": [
    "# PGA"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "a61720f7",
   "metadata": {},
   "outputs": [],
   "source": [
    "PGA is not an independent assembly software, it relies on the output of hifiasm (–write-ec –write-paf -l 1) and technologies that can anchor contigs to chromosomes,such as genetic maps, HIC, or alignment with genome sequences of related species to achieve a gapless genome."
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "f7125fc8",
   "metadata": {},
   "outputs": [],
   "source": [
    "# Dependencies"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "aa8e56a0",
   "metadata": {},
   "outputs": [],
   "source": [
    "python3\n",
    "\n",
    "mummer\n",
    "\n",
    "networkx"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "67dc975c",
   "metadata": {},
   "outputs": [],
   "source": [
    "# Installation"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "a9f9f596",
   "metadata": {},
   "outputs": [],
   "source": [
    "conda install -c bioconda mummer\n",
    "\n",
    "conda install -c conda-forge networkx\n",
    "\n",
    "git clone https://github.com/likui345/PGA.git\n",
    "\n",
    "python setup.py install"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "a728beda",
   "metadata": {},
   "outputs": [],
   "source": [
    "# Usage"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "072e9d84",
   "metadata": {},
   "outputs": [],
   "source": [
    "### 1.Using reference genomes to anchor scaffolds."
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "1232b9ce",
   "metadata": {},
   "outputs": [],
   "source": [
    "\n",
    "usage: \n",
    "    rbsa.py  --type nucmer --ref REF --scf SCF\n",
    "    or rbsa.py  --type mcscan --A_bed A_BED --B_bed B_BED --anchor ANCHOR --scf SCF\n",
    "    \n",
    "optional arguments:\n",
    "  -h, --help            show this help message and exit\n",
    "  --type {nucmer,mcscan}\n",
    "                        nucmer or mcscan anchor\n",
    "  --ref REF             This is reference genome sequences which has anchored\n",
    "                        to the chromosomal level\n",
    "  --scf SCF             This is the scaffold sequences\n",
    "  --A_bed A_BED         species_a bed file\n",
    "  --B_bed B_BED         species_b bed file\n",
    "  --anchor ANCHOR       anchor.simple files of species_a and species_b\n",
    "    "
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "60551edb",
   "metadata": {},
   "outputs": [],
   "source": [
    "### 2.Filtering overlap data ."
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "634d5c04",
   "metadata": {},
   "outputs": [],
   "source": [
    "\n",
    "usage: \n",
    "    filter.py  --paf_fn PAF --bestn NUM --output OUTPUT\n",
    "    \n",
    "optional arguments:\n",
    "  -h, --help       show this help message and exit\n",
    "  --paf_fn PAF_FN  Input. reads alignment file\n",
    "  --bestn BESTN    output at least best n overlaps on 5' or 3' ends if\n",
    "                   possible.\n",
    "  --output OUTPUT  Output filename"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "2e32d996",
   "metadata": {},
   "outputs": [],
   "source": [
    "### 3.Get chr_paths using agp file and gfa file."
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "3272a750",
   "metadata": {},
   "outputs": [],
   "source": [
    "\n",
    "usage: chr_paths.py  --agp AGP --gfa GFA\n",
    "\n",
    "optional arguments:\n",
    "  -h, --help  show this help message and exit\n",
    "  --agp AGP   agp file form contig anchoring\n",
    "  --gfa GFA   gfa from hifiasm"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "171cef32",
   "metadata": {},
   "outputs": [],
   "source": [
    "### 4.Building string graph."
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "6eeff565",
   "metadata": {},
   "outputs": [],
   "source": [
    "\n",
    "usage: ovlp2graph.py [-h] --overlap-file OVERLAP_FILE\n",
    "\n",
    "string graph assembler that is desinged for handling diploid genomes\n",
    "\n",
    "optional arguments:\n",
    "  -h, --help            show this help message and exit\n",
    "  --overlap-file OVERLAP_FILE\n",
    "                        the filtered overlap data from step2. (default: None)\n",
    "\n",
    "Outputs:\n",
    "    - sg_edges_list\n",
    "    - chimer_nodes \n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "abebc7c1",
   "metadata": {},
   "outputs": [],
   "source": [
    "### 5.Generate the final assembly."
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "7090915c",
   "metadata": {},
   "outputs": [],
   "source": [
    "\n",
    "usage: graph2chr.py [-h] --reads-fasta-fn READS_FASTA_FN --paf-fn PAF_FN\n",
    "                     --sg-edges-list-fn SG_EDGES_LIST_FN --chr-paths-fn\n",
    "                     CHR_PATHS_FN\n",
    "\n",
    "Generate the chromosome and tiling paths, given the string graph.\n",
    "\n",
    "optional arguments:\n",
    "  -h, --help            show this help message and exit\n",
    "  --reads-fasta-fn READS_FASTA_FN\n",
    "                        Input. reads file generated by hifiasm. (default:\n",
    "                        None)\n",
    "  --paf-fn PAF_FN       Input. reads alignment file generated by hifiasm.\n",
    "                        (default: None)\n",
    "  --sg-edges-list-fn SG_EDGES_LIST_FN\n",
    "                        Input. File containing string graph edges, produced by\n",
    "                        ovlp2graph.py. (default: None)\n",
    "  --chr-paths-fn CHR_PATHS_FN\n",
    "                        Input. File containing chromosome paths. (default:\n",
    "                        None)\n",
    "\n",
    "\n",
    "\n"
   ]
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": "Python 3",
   "language": "python",
   "name": "python3"
  },
  "language_info": {
   "codemirror_mode": {
    "name": "ipython",
    "version": 3
   },
   "file_extension": ".py",
   "mimetype": "text/x-python",
   "name": "python",
   "nbconvert_exporter": "python",
   "pygments_lexer": "ipython3",
   "version": "3.8.8"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 5
}
