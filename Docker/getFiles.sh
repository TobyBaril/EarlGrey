#!/bin/sh

set -eu

# download src name
download() {
  src="$1"
  shift

  if [ $# -ge 1 ]; then
    name="$1"
  else
    name="${src##*/}"
  fi

  dest="src/$name"

  if [ -n "${ALWAYS-}" ] || ! [ -f "$dest" ]; then
    echo "Downloading $src to $dest"
    curl -sSL "$src" > "$dest"
  fi
}

mkdir -p src

download https://www.repeatmasker.org/rmblast/rmblast-2.13.0+-x64-linux.tar.gz
download http://eddylab.org/software/hmmer/hmmer-3.3.2.tar.gz
download https://github.com/Benson-Genomics-Lab/TRF/archive/v4.09.1.tar.gz trf-4.09.1.tar.gz
download https://www.repeatmasker.org/RepeatScout-1.0.6.tar.gz
download https://www.repeatmasker.org/RepeatModeler/RECON-1.08.tar.gz
download https://github.com/weizhongli/cdhit/releases/download/V4.8.1/cd-hit-v4.8.1-2019-0228.tar.gz
download https://github.com/genometools/genometools/archive/v1.6.2.tar.gz gt-1.6.2.tar.gz
download https://github.com/oushujun/LTR_retriever/archive/v2.9.0.tar.gz LTR_retriever-2.9.0.tar.gz
download https://mafft.cbrc.jp/alignment/software/mafft-7.471-without-extensions-src.tgz
download https://github.com/TravisWheelerLab/NINJA/archive/0.97-cluster_only.tar.gz NINJA-cluster.tar.gz
download https://www.repeatmasker.org/coseg-0.2.2.tar.gz
download https://www.repeatmasker.org/RepeatMasker/RepeatMasker-4.1.4.tar.gz
download https://github.com/Dfam-consortium/RepeatModeler/archive/2.0.4.tar.gz RepeatModeler-2.0.4.tar.gz
download https://github.com/Dfam-consortium/FamDB/raw/master/famdb.py famdb.py
download https://www.dfam.org/releases/current/families/Dfam_curatedonly.h5.gz Dfam.h5.gz

# TODO: /exe/ only includes binaries of the "latest" version at the time of download.
# The version listed in README.md is obtained by running 'strings src/faToTwoBit | grep kent'
# On whatever was downloaded.
# Consider building these tools from source instead.
for tool in faToTwoBit twoBitInfo twoBitToFa; do
  download https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/"$tool"
done
