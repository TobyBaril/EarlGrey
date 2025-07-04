#!/bin/bash

usage() { echo "Usage: [-l Repeat library] [-g Genome ] [-t Threads (default 4) ] [-f Flank (default 1000) ] [-m Minimum number of sequences required in multiple sequence alignments during BEAT] [-r Number of iterations of BEAT to run (deafult 10)] [-d Out directory, if not specified wil be created ] [-h Print this help] [-M Amount of memory TEstrainer needs to keep free] [-R Resume interrupted run in directory specified by -d]" 1>&2; exit 1; }

# default values
STRAIN_SCRIPTS=/data/toby/EarlGrey/scripts/TEstrainer/scripts/
FLANK=1000
THREADS=4
RUNS=10
NO_SEQ=20
# for potential folder name
TIME=$(date +"%s")
TIME=${TIME: -4}
MEM_FREE="200M"
MIN_SEQ=3

# parsing
while getopts l:g:t:f:m:r:d:h:n:M:R flag; do
  case "${flag}" in
    l) RM_LIBRARY_PATH=${OPTARG};;
    g) GENOME=${OPTARG};;
    t) THREADS=${OPTARG};;
    f) FLANK=${OPTARG};;
    m) MIN_SEQ=${OPTARG};;
    r) RUNS=${OPTARG};;
    d) DATA_DIR=${OPTARG};;
    n) NO_SEQ=${OPTARG};;
    M) MEM_FREE=${OPTARG};;
    R) RESUME=true;;
    h | *)
      print_usage
      exit_script
  esac
done

# determine further variables and check files exist
if [ ! -f ${RM_LIBRARY_PATH} ]; then echo "Library not found"; usage; fi
if [ -z ${RM_LIBRARY_PATH} ]; then echo "Library must be supplied"; usage; else RM_LIBRARY=$(echo $RM_LIBRARY_PATH | sed 's/.*\///'); fi
if [[ $RUNS -gt 0 ]]; then 
  if [ -z ${GENOME} ]; then echo "If refining genome must be supplied"; usage; fi
  if [ ! -f ${GENOME} ]; then echo "Refining genome not found"; usage; fi
fi

if [[ $THREADS -gt 4 ]]; then MAFFT_THREADS=$(($(($THREADS / 4)))); else MAFFT_THREADS=1; fi

if [ -n "$RESUME" ]
then
  if [ -z "$DATA_DIR" ]; then echo "When resuming, out directory (-d) must be specified"; usage; fi
  RUN_NO=$(ls ${DATA_DIR} | grep -P '^run_\d+$' | sed 's/run_//' | sort -n | tail -n 1)
  if [ -z "$RUN_NO" ] || [ $RUN_NO -eq 0 ]
  then
    echo "Attempting to resume a run that was interrupted before iterations started"
    unset RESUME
    if [ -d ${DATA_DIR}/run_0/ ]; then rm -r ${DATA_DIR}/run_0/ ; fi
    if [ -d ${DATA_DIR}/tmp/ ]; then rm -r ${DATA_DIR}/tmp/ ; fi
  elif [ $RUN_NO -gt $RUNS ]
  then
    echo "Trying to resume, but more runs than expected have already been completed"
    usage
  elif [ $RUN_NO -eq $RUNS ] && [ -f ${DATA_DIR}/${RM_LIBRARY} ]
  then
    RUN_NO=$RUN_NO+1
  else
    rm -r ${DATA_DIR}/run_${RUN_NO}
  fi
fi

if [ -z "$RESUME" ]
then
  if [ -z "$DATA_DIR" ]; then DATA_DIR=$(echo "TS_"${RM_LIBRARY}"_"${TIME}); fi

  # create data dir if missing
  if [ ! -d ${DATA_DIR} ]
  then
    mkdir ${DATA_DIR}
  fi

  # make directories
  mkdir -p ${DATA_DIR}/run_0/
  mkdir -p ${DATA_DIR}/tmp/

  # initial copy
  cp ${RM_LIBRARY_PATH} ${DATA_DIR}/${RM_LIBRARY}

  # cp starting seq to starting directory
  mkdir -p ${DATA_DIR}/run_0/og
  cp ${RM_LIBRARY_PATH} ${DATA_DIR}/run_0/further_${RM_LIBRARY}
  # create reference of original sequences
  python ${STRAIN_SCRIPTS}/splitter.py -i ${DATA_DIR}/run_0/further_${RM_LIBRARY} -o ${DATA_DIR}/run_0/og

  # runs
  python ${STRAIN_SCRIPTS}/indexer.py -g ${GENOME}
  if [ ! -f "${GENOME}".nsq ]; then
    makeblastdb -in ${GENOME} -dbtype nucl -out ${GENOME} # makeblastb if needed
  fi

  # curation
  RUN_NO=1
fi

while  [ $RUN_NO -le $RUNS ]
do

  # make directories
  mkdir -p ${DATA_DIR}/run_${RUN_NO}/raw \
           ${DATA_DIR}/run_${RUN_NO}/TEtrim_complete \
           ${DATA_DIR}/run_${RUN_NO}/initial_blast \
           ${DATA_DIR}/run_${RUN_NO}/self_search \
           ${DATA_DIR}/run_${RUN_NO}/to_align \
           ${DATA_DIR}/run_${RUN_NO}/mafft \
           ${DATA_DIR}/run_${RUN_NO}/TEtrim_con \
           ${DATA_DIR}/run_${RUN_NO}/TEtrim_unaln \
           ${DATA_DIR}/run_${RUN_NO}/TEtrim_blast \
           ${DATA_DIR}/run_${RUN_NO}/TEtrim_mafft \
           ${DATA_DIR}/run_${RUN_NO}/TEtrim_further \
           ${DATA_DIR}/run_${RUN_NO}/TEtrim_bp
  
  # split
  cp ${DATA_DIR}/run_$(expr $RUN_NO - 1)/further_${RM_LIBRARY} ${DATA_DIR}/run_${RUN_NO}/${RM_LIBRARY}
  echo "Splitting run "${RUN_NO}
  python ${STRAIN_SCRIPTS}/splitter.py -i ${DATA_DIR}/run_${RUN_NO}/${RM_LIBRARY} -o ${DATA_DIR}/run_${RUN_NO}/raw

  # run trf to determine if sequence is tandem repeat
  echo "Initial trf check for "${RUN_NO}
  parallel --bar --tmpdir ${DATA_DIR}/tmp/ --jobs ${THREADS} --memfree ${MEM_FREE} -a ${DATA_DIR}/run_${RUN_NO}/raw/${RM_LIBRARY}_split.txt trf ${DATA_DIR}/run_${RUN_NO}/raw/{} 2 7 7 80 10 50 500 -d -h -ngs ">" ${DATA_DIR}/run_${RUN_NO}/raw/{}.trf
  echo "Initial blast and preparation for MSA "${RUN_NO}
  # initial blast and extention
  parallel --bar --tmpdir ${DATA_DIR}/tmp/ --jobs ${THREADS} --memfree ${MEM_FREE} -a ${DATA_DIR}/run_${RUN_NO}/raw/${RM_LIBRARY}_split.txt python3 ${STRAIN_SCRIPTS}/initial_mafft_setup.py -d ${DATA_DIR} -r ${RUN_NO} -s {} -g ${GENOME} -f ${FLANK} -D -n ${NO_SEQ}
  
  ## first mafft alignment
  find ${DATA_DIR}/run_${RUN_NO}/to_align -type f | sed 's/.*\///' > ${DATA_DIR}/run_${RUN_NO}/to_align.txt
  echo "Primary alignment run "${RUN_NO}
  parallel --bar --tmpdir ${DATA_DIR}/tmp/ --jobs $MAFFT_THREADS -a ${DATA_DIR}/run_${RUN_NO}/to_align.txt mafft --thread 4 --quiet --localpair --adjustdirectionaccurately ${DATA_DIR}/run_${RUN_NO}/to_align/{} ">" ${DATA_DIR}/run_${RUN_NO}/mafft/{}

  # trim
  echo "Trimming run "${RUN_NO}
  parallel --bar --tmpdir ${DATA_DIR}/tmp/ --jobs $MAFFT_THREADS -a ${DATA_DIR}/run_${RUN_NO}/to_align.txt python ${STRAIN_SCRIPTS}/TEtrim.py -i ${DATA_DIR}/run_${RUN_NO}/mafft/{} -t 4 -f ${FLANK} -n ${RUN_NO} -d ${DATA_DIR} -m ${MIN_SEQ}
  
  # compile completed curations
  if [ -n "$(ls -A ${DATA_DIR}/run_${RUN_NO}/TEtrim_complete/ 2>/dev/null)" ]; then
     cat ${DATA_DIR}/run_${RUN_NO}/TEtrim_complete/*fasta > ${DATA_DIR}/run_${RUN_NO}/complete_${RM_LIBRARY}
  fi
  
  # Either compile consensuses for further curation
  if [ -n "$(ls -A ${DATA_DIR}/run_${RUN_NO}/TEtrim_further/ 2>/dev/null)" ]; then
    cat ${DATA_DIR}/run_${RUN_NO}/TEtrim_further/*fasta > ${DATA_DIR}/run_${RUN_NO}/further_${RM_LIBRARY}
    echo "Ready for " $(( RUN_NO++ ))
  else
  # or exit loop if no more curation needed
    echo "Finished extension"
    break
  fi

done

# Compile all completed
cat ${DATA_DIR}/run_*/complete_${RM_LIBRARY} > ${DATA_DIR}/${RM_LIBRARY}
# if more could have been extended append to compilation
if [ -s ${DATA_DIR}/run_${RUNS}/further_${RM_LIBRARY} ]; then
  cat ${DATA_DIR}/run_${RUNS}/further_${RM_LIBRARY} >> ${DATA_DIR}/${RM_LIBRARY}
fi

# add any families which went missing along the way due to alignment issues
find ${DATA_DIR}/run_*/mafft/ -empty -type f | sed 's/mafft/raw/' > ${DATA_DIR}/missing_consensi.txt
if [[ -s ${DATA_DIR}/missing_consensi.txt ]]; then
  while read file_name; do
    cat $file_name >> ${DATA_DIR}/${RM_LIBRARY}
  done < ${DATA_DIR}/missing_consensi.txt
fi
  
sed -i 's/ .*//' ${DATA_DIR}/${RM_LIBRARY}

# Identify simple repeats and satellites, trim ends of LINEs/SINEs
echo "Splitting for simple/satellite packages"
mkdir -p ${DATA_DIR}/trf/split
python ${STRAIN_SCRIPTS}/splitter.py -i ${DATA_DIR}/${RM_LIBRARY} -o ${DATA_DIR}/trf/split
cp ${DATA_DIR}/${RM_LIBRARY} ${DATA_DIR}/trf/
# Run and parse TRF
echo "Running TRF"
parallel --bar --tmpdir ${DATA_DIR}/tmp/ --jobs ${THREADS} --memfree ${MEM_FREE} -a ${DATA_DIR}/trf/split/${RM_LIBRARY}_split.txt trf ${DATA_DIR}/trf/split/{} 2 7 7 80 10 50 500 -d -h -ngs ">" ${DATA_DIR}/trf/split/{}.trf
parallel --bar --tmpdir ${DATA_DIR}/tmp/ --jobs ${THREADS} --memfree ${MEM_FREE} -a ${DATA_DIR}/trf/split/${RM_LIBRARY}_split.txt python3 ${STRAIN_SCRIPTS}/trf_parser.py --trf ${DATA_DIR}/trf/split/{}.trf --out ${DATA_DIR}/trf/split/{}.trf.tsv
find ${DATA_DIR}/trf/split/ -type f -name "*trf.tsv" -exec cat {} + | cat > ${DATA_DIR}/trf/${RM_LIBRARY}.trf
# Run SA-SSR
echo "Running SA-SSR"
sa-ssr -e -l 20 -L 50000 -m 1 -M 5000 -t ${THREADS} ${DATA_DIR}/${RM_LIBRARY} ${DATA_DIR}/trf/${RM_LIBRARY}.sassr
# Run and compile mreps
echo "Running mreps"
parallel --bar --tmpdir ${DATA_DIR}/tmp/ --jobs ${THREADS} --memfree ${MEM_FREE} -a ${DATA_DIR}/trf/split/${RM_LIBRARY}_split.txt bash ${STRAIN_SCRIPTS}/mreps_parser.sh -i ${DATA_DIR}/trf/split/{} "2>"/dev/null
find ${DATA_DIR}/trf/split/ -type f -name "*mreps" -exec cat {} + | cat > ${DATA_DIR}/trf/${RM_LIBRARY}.mreps
# Interpret mreps, TRF and SA-SSR
echo "Trimming and sorting based on mreps, TRF, SA-SSR"
if [ ! -f ${DATA_DIR}/trf/${RM_LIBRARY}.sassr ]; then
   touch ${DATA_DIR}/trf/${RM_LIBRARY}.sassr
fi
Rscript ${STRAIN_SCRIPTS}/simple_repeat_filter_trim.R -i ${DATA_DIR}/${RM_LIBRARY} -d ${DATA_DIR} -p 50

# Delete temp files
echo "Removing temporary files"
rm -r ${DATA_DIR}/*/split/
rm -r ${DATA_DIR}/tmp/
find ${DATA_DIR}/ -mindepth 1 -name "run_*" -exec rm -r {} +
if [[ $RUNS -gt 0 ]]; then rm ${GENOME}.n*; fi

# Classify improved consensi using RepeatModeler's RepeatClassifier
echo "Reclassifying repeats"
mkdir -p ${DATA_DIR}/classify/
cp ${DATA_DIR}/trf/${RM_LIBRARY}.nonsatellite ${DATA_DIR}/classify/
cd ${DATA_DIR}/classify/
if [ -s ${RM_LIBRARY}.nonsatellite ]; then
    RepeatClassifier -threads ${THREADS} -consensi ${RM_LIBRARY}.nonsatellite
fi
# legacy command for older installations - DEPRECATED
# RepeatClassifier -pa ${THREADS} -consensi ${RM_LIBRARY}.nonsatellite
cd -
# Compile classified files
if [ -f ${DATA_DIR}/classify/${RM_LIBRARY}.nonsatellite.classified ]; then
    cp ${DATA_DIR}/classify/${RM_LIBRARY}.nonsatellite.classified ${DATA_DIR}/${RM_LIBRARY}.strained
else
    touch ${DATA_DIR}/${RM_LIBRARY}.strained
fi
echo "Compiling library"
if [ -f ${DATA_DIR}/trf/${RM_LIBRARY}.satellites ]; then
    cat ${DATA_DIR}/trf/${RM_LIBRARY}.satellites >> ${DATA_DIR}/${RM_LIBRARY}.strained
fi
