#!/bin/bash
set -eo pipefail

usage() { echo "Usage: [-l Repeat library] [-g Genome] [-t Threads (default 4)] [-f Flank (default 1000)] [-m Minimum sequences] [-r Iterations (default 10)] [-d Output directory] [-n Maximum hits] [-B Memory budget in GiB (default 64)] [-J Maximum concurrent MAFFT jobs (default auto)] [-A Alignment mode: localpair or fast (default localpair)] [-M Memory GNU parallel should keep free] [-R Resume] [-q Suppress progress bar] [-h Help]" 1>&2; exit 1; }

# default values
STRAIN_SCRIPTS=/data/toby/EarlGrey/scripts/TEstrainer/scripts/
FLANK=1000
THREADS=4
RUNS=10
NO_SEQ=20
# for potential folder name
TIME=$(date +"%s")
TIME=${TIME: -4}
MEM_FREE="1G"
MIN_SEQ=3
BAR_FLAG="--bar"
MEMORY_GIB=64
ALIGN_JOBS=auto
ALIGN_MODE=localpair

# parsing
while getopts l:g:t:f:m:r:d:h:n:M:B:J:A:Rq flag; do
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
    B) MEMORY_GIB=${OPTARG};;
    J) ALIGN_JOBS=${OPTARG};;
    A) ALIGN_MODE=${OPTARG};;
    R) RESUME=true;;
    q) BAR_FLAG="";;
    h | *)
      usage
  esac
done
case "$ALIGN_MODE" in
  localpair) PRIMARY_FLAGS=(--large --localpair --adjustdirectionaccurately); TRIM_FLAGS=(--large);;
  fast) PRIMARY_FLAGS=(--memsave --6merpair --adjustdirection); TRIM_FLAGS=(--fast);;
  *) echo "ERROR: alignment mode must be localpair or fast" >&2; exit 1;;
esac

# determine further variables and check files exist
if [ ! -f ${RM_LIBRARY_PATH} ]; then echo "Library not found"; usage; fi
if [ -z ${RM_LIBRARY_PATH} ]; then echo "Library must be supplied"; usage; else RM_LIBRARY=$(echo $RM_LIBRARY_PATH | sed 's/.*\///'); fi
if [[ $RUNS -gt 0 ]]; then 
  if [ -z ${GENOME} ]; then echo "If refining genome must be supplied"; usage; fi
  if [ ! -f ${GENOME} ]; then echo "Refining genome not found"; usage; fi
fi

for value in "$THREADS" "$MEMORY_GIB"; do
  if ! [[ "$value" =~ ^[1-9][0-9]*$ ]]; then
    echo "ERROR: threads and memory budget must be positive integers" >&2
    exit 1
  fi
done
if [[ "$ALIGN_JOBS" != auto ]] && ! [[ "$ALIGN_JOBS" =~ ^[1-9][0-9]*$ ]]; then
  echo "ERROR: alignment jobs must be a positive integer or auto" >&2
  exit 1
fi
GUARD=${STRAIN_SCRIPTS}/memory_guard.py
USABLE_MEM_MB=$(python3 "$GUARD" --budget "$((MEMORY_GIB * 1024))")
MAX_JOBS=$((USABLE_MEM_MB / 1024))
if (( MAX_JOBS < 1 )); then
  echo "ERROR: less than 1 GiB is available for TEstrainer subprocesses" >&2
  exit 1
fi
PREP_JOBS=$THREADS
if (( PREP_JOBS > 4 )); then PREP_JOBS=4; fi
if (( PREP_JOBS > MAX_JOBS )); then PREP_JOBS=$MAX_JOBS; fi
if [[ "$ALIGN_JOBS" == auto ]]; then
  ALIGN_JOBS=$((USABLE_MEM_MB / 8192))
  if (( ALIGN_JOBS < 1 )); then ALIGN_JOBS=1; fi
fi
if (( ALIGN_JOBS > THREADS )); then ALIGN_JOBS=$THREADS; fi
if (( ALIGN_JOBS > MAX_JOBS )); then ALIGN_JOBS=$MAX_JOBS; fi
PREP_MEM_MB=$((USABLE_MEM_MB / PREP_JOBS))
ALIGN_MEM_MB=$((USABLE_MEM_MB / ALIGN_JOBS))
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
echo "TEstrainer memory: ${USABLE_MEM_MB} MiB usable; ${PREP_JOBS} preparation jobs (${PREP_MEM_MB} MiB each), ${ALIGN_JOBS} alignment jobs (${ALIGN_MEM_MB} MiB each; 1 MAFFT thread each; mode ${ALIGN_MODE})"

if [ -n "$RESUME" ]
then
  if [ -z "$DATA_DIR" ]; then echo "When resuming, out directory (-d) must be specified"; usage; fi
  if [ ! -d "$DATA_DIR" ]; then echo "ERROR: resume directory not found: $DATA_DIR" >&2; exit 1; fi
  PREVIOUS_MODE=localpair
  if [ -f "${DATA_DIR}/alignment_mode" ]; then
    PREVIOUS_MODE=$(cat "${DATA_DIR}/alignment_mode")
  fi
  if [ "$PREVIOUS_MODE" != "$ALIGN_MODE" ]; then
    echo "ERROR: cannot resume a ${PREVIOUS_MODE} run using ${ALIGN_MODE}; start a fresh run" >&2
    exit 1
  fi
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
  elif [ -f "${DATA_DIR}/run_${RUN_NO}/.complete" ]
  then
    if [ -s "${DATA_DIR}/run_${RUN_NO}/further_${RM_LIBRARY}" ]; then
      RUN_NO=$((RUN_NO + 1))
    else
      RUN_NO=$((RUNS + 1))
    fi
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
  python3 ${GUARD} --limit ${USABLE_MEM_MB} --label initial-split -- python ${STRAIN_SCRIPTS}/splitter.py -i ${DATA_DIR}/run_0/further_${RM_LIBRARY} -o ${DATA_DIR}/run_0/og

  # runs
  python3 ${GUARD} --limit ${USABLE_MEM_MB} --label genome-index -- python ${STRAIN_SCRIPTS}/indexer.py -g ${GENOME}
  if [ ! -f "${GENOME}".nsq ]; then
    python3 ${GUARD} --limit ${USABLE_MEM_MB} --label makeblastdb -- makeblastdb -in ${GENOME} -dbtype nucl -out ${GENOME}
  fi

  # curation
  RUN_NO=1
fi

printf '%s\n' "$ALIGN_MODE" > "${DATA_DIR}/alignment_mode"
mkdir -p "${DATA_DIR}/tmp/"
if [ -z "${MAFFT_TMPDIR:-}" ]; then
  MAFFT_TMPDIR=$(cd "${DATA_DIR}/tmp/" && pwd)
fi
mkdir -p "$MAFFT_TMPDIR"
export MAFFT_TMPDIR
echo "MAFFT scratch directory: $MAFFT_TMPDIR"

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
  python3 ${GUARD} --limit ${USABLE_MEM_MB} --label split-run-${RUN_NO} -- python ${STRAIN_SCRIPTS}/splitter.py -i ${DATA_DIR}/run_${RUN_NO}/${RM_LIBRARY} -o ${DATA_DIR}/run_${RUN_NO}/raw

  # run trf to determine if sequence is tandem repeat
  echo "Initial trf check for "${RUN_NO}
  parallel ${BAR_FLAG} --halt soon,fail=1 --tmpdir ${DATA_DIR}/tmp/ --jobs ${PREP_JOBS} --memfree ${MEM_FREE} -a ${DATA_DIR}/run_${RUN_NO}/raw/${RM_LIBRARY}_split.txt python3 ${GUARD} --limit ${PREP_MEM_MB} --label trf:{} --output ${DATA_DIR}/run_${RUN_NO}/raw/{}.trf --allow-empty-output --fallback-source ${DATA_DIR}/run_${RUN_NO}/raw/{} --fallback-output ${DATA_DIR}/run_${RUN_NO}/TEtrim_complete/{} -- trf ${DATA_DIR}/run_${RUN_NO}/raw/{} 2 7 7 80 10 50 500 -d -h -ngs
  echo "Initial blast and preparation for MSA "${RUN_NO}
  # initial blast and extention
  while IFS= read -r family; do
    if [ ! -f "${DATA_DIR}/run_${RUN_NO}/TEtrim_complete/${family}" ]; then printf '%s\n' "$family"; fi
  done < "${DATA_DIR}/run_${RUN_NO}/raw/${RM_LIBRARY}_split.txt" > "${DATA_DIR}/run_${RUN_NO}/to_prepare.txt"
  parallel ${BAR_FLAG} --halt soon,fail=1 --tmpdir ${DATA_DIR}/tmp/ --jobs ${PREP_JOBS} --memfree ${MEM_FREE} -a ${DATA_DIR}/run_${RUN_NO}/to_prepare.txt python3 ${GUARD} --limit ${PREP_MEM_MB} --label preparation:{} --fallback-source ${DATA_DIR}/run_${RUN_NO}/raw/{} --fallback-output ${DATA_DIR}/run_${RUN_NO}/TEtrim_complete/{} --remove-on-fallback ${DATA_DIR}/run_${RUN_NO}/to_align/{} -- python3 ${STRAIN_SCRIPTS}/initial_mafft_setup.py -d ${DATA_DIR} -r ${RUN_NO} -s {} -g ${GENOME} -f ${FLANK} -D -n ${NO_SEQ}
  
  ## first mafft alignment
  find ${DATA_DIR}/run_${RUN_NO}/to_align -type f | sed 's/.*\///' > ${DATA_DIR}/run_${RUN_NO}/to_align.txt
  echo "Primary alignment run "${RUN_NO}
  parallel ${BAR_FLAG} --halt soon,fail=1 --tmpdir ${DATA_DIR}/tmp/ --jobs ${ALIGN_JOBS} --memfree ${MEM_FREE} -a ${DATA_DIR}/run_${RUN_NO}/to_align.txt python3 ${GUARD} --limit ${ALIGN_MEM_MB} --label primary-mafft:{} --output ${DATA_DIR}/run_${RUN_NO}/mafft/{} --fallback-source ${DATA_DIR}/run_${RUN_NO}/raw/{} --fallback-output ${DATA_DIR}/run_${RUN_NO}/TEtrim_complete/{} --remove-on-fallback ${DATA_DIR}/run_${RUN_NO}/mafft/{} -- mafft --thread 1 --quiet "${PRIMARY_FLAGS[@]}" ${DATA_DIR}/run_${RUN_NO}/to_align/{}

  # trim
  echo "Trimming run "${RUN_NO}
  while IFS= read -r family; do
    if [ -f "${DATA_DIR}/run_${RUN_NO}/mafft/${family}" ]; then printf '%s\n' "$family"; fi
  done < "${DATA_DIR}/run_${RUN_NO}/to_align.txt" > "${DATA_DIR}/run_${RUN_NO}/to_trim.txt"
  parallel ${BAR_FLAG} --halt soon,fail=1 --tmpdir ${DATA_DIR}/tmp/ --jobs ${ALIGN_JOBS} --memfree ${MEM_FREE} -a ${DATA_DIR}/run_${RUN_NO}/to_trim.txt python3 ${GUARD} --limit ${ALIGN_MEM_MB} --label trimming:{} --fallback-source ${DATA_DIR}/run_${RUN_NO}/raw/{} --fallback-output ${DATA_DIR}/run_${RUN_NO}/TEtrim_complete/{} --remove-on-fallback ${DATA_DIR}/run_${RUN_NO}/TEtrim_further/{} -- python ${STRAIN_SCRIPTS}/TEtrim.py -i ${DATA_DIR}/run_${RUN_NO}/mafft/{} -t 1 -f ${FLANK} -n ${RUN_NO} -d ${DATA_DIR} -m ${MIN_SEQ} "${TRIM_FLAGS[@]}"
  
  # compile completed curations
  if [ -n "$(ls -A ${DATA_DIR}/run_${RUN_NO}/TEtrim_complete/ 2>/dev/null)" ]; then
     cat ${DATA_DIR}/run_${RUN_NO}/TEtrim_complete/*fasta > ${DATA_DIR}/run_${RUN_NO}/complete_${RM_LIBRARY}
  fi
  
  # Either compile consensuses for further curation
  if [ -n "$(ls -A ${DATA_DIR}/run_${RUN_NO}/TEtrim_further/ 2>/dev/null)" ]; then
    cat ${DATA_DIR}/run_${RUN_NO}/TEtrim_further/*fasta > ${DATA_DIR}/run_${RUN_NO}/further_${RM_LIBRARY}
    touch "${DATA_DIR}/run_${RUN_NO}/.complete"
    echo "Ready for " "${RUN_NO}"
    RUN_NO=$(( RUN_NO + 1 ))
  else
  # or exit loop if no more curation needed
    touch "${DATA_DIR}/run_${RUN_NO}/.complete"
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
python3 ${GUARD} --limit ${USABLE_MEM_MB} --label final-split -- python ${STRAIN_SCRIPTS}/splitter.py -i ${DATA_DIR}/${RM_LIBRARY} -o ${DATA_DIR}/trf/split
cp ${DATA_DIR}/${RM_LIBRARY} ${DATA_DIR}/trf/
# Run and parse TRF
echo "Running TRF"
parallel ${BAR_FLAG} --halt soon,fail=1 --tmpdir ${DATA_DIR}/tmp/ --jobs ${PREP_JOBS} --memfree ${MEM_FREE} -a ${DATA_DIR}/trf/split/${RM_LIBRARY}_split.txt python3 ${GUARD} --limit ${PREP_MEM_MB} --label trf:{} -- trf ${DATA_DIR}/trf/split/{} 2 7 7 80 10 50 500 -d -h -ngs ">" ${DATA_DIR}/trf/split/{}.trf
parallel ${BAR_FLAG} --halt soon,fail=1 --tmpdir ${DATA_DIR}/tmp/ --jobs ${PREP_JOBS} --memfree ${MEM_FREE} -a ${DATA_DIR}/trf/split/${RM_LIBRARY}_split.txt python3 ${GUARD} --limit ${PREP_MEM_MB} --label trf-parser:{} -- python3 ${STRAIN_SCRIPTS}/trf_parser.py --trf ${DATA_DIR}/trf/split/{}.trf --out ${DATA_DIR}/trf/split/{}.trf.tsv
find ${DATA_DIR}/trf/split/ -type f -name "*trf.tsv" -exec cat {} + | cat > ${DATA_DIR}/trf/${RM_LIBRARY}.trf
# Run SA-SSR
echo "Running SA-SSR"
python3 ${GUARD} --limit ${USABLE_MEM_MB} --label sa-ssr -- sa-ssr -e -l 20 -L 50000 -m 1 -M 5000 -t ${THREADS} ${DATA_DIR}/${RM_LIBRARY} ${DATA_DIR}/trf/${RM_LIBRARY}.sassr
# Run and compile mreps
echo "Running mreps"
parallel ${BAR_FLAG} --halt soon,fail=1 --tmpdir ${DATA_DIR}/tmp/ --jobs ${PREP_JOBS} --memfree ${MEM_FREE} -a ${DATA_DIR}/trf/split/${RM_LIBRARY}_split.txt python3 ${GUARD} --limit ${PREP_MEM_MB} --label mreps:{} -- bash ${STRAIN_SCRIPTS}/mreps_parser.sh -i ${DATA_DIR}/trf/split/{}
find ${DATA_DIR}/trf/split/ -type f -name "*mreps" -exec cat {} + | cat > ${DATA_DIR}/trf/${RM_LIBRARY}.mreps
# Interpret mreps, TRF and SA-SSR
echo "Trimming and sorting based on mreps, TRF, SA-SSR"
if [ ! -f ${DATA_DIR}/trf/${RM_LIBRARY}.sassr ]; then
   touch ${DATA_DIR}/trf/${RM_LIBRARY}.sassr
fi
python3 ${GUARD} --limit ${USABLE_MEM_MB} --label repeat-filter -- Rscript ${STRAIN_SCRIPTS}/simple_repeat_filter_trim.R -i ${DATA_DIR}/${RM_LIBRARY} -d ${DATA_DIR} -p 50

# Classify improved consensi using RepeatModeler's RepeatClassifier
echo "Reclassifying repeats"
mkdir -p ${DATA_DIR}/classify/
cp ${DATA_DIR}/trf/${RM_LIBRARY}.nonsatellite ${DATA_DIR}/classify/
cd ${DATA_DIR}/classify/
if [ -s ${RM_LIBRARY}.nonsatellite ]; then
    rm -f "${RM_LIBRARY}.nonsatellite.classified"
    python3 ${GUARD} --limit ${USABLE_MEM_MB} --label RepeatClassifier -- RepeatClassifier -threads ${THREADS} -consensi ${RM_LIBRARY}.nonsatellite
    if [ ! -s ${RM_LIBRARY}.nonsatellite.classified ]; then
      echo "ERROR: RepeatClassifier did not produce a classified library" >&2
      exit 1
    fi
fi
# legacy command for older installations - DEPRECATED
# RepeatClassifier -pa ${THREADS} -consensi ${RM_LIBRARY}.nonsatellite
cd -
# Compile classified files
if [ -s ${DATA_DIR}/trf/${RM_LIBRARY}.nonsatellite ]; then
    cp ${DATA_DIR}/classify/${RM_LIBRARY}.nonsatellite.classified ${DATA_DIR}/${RM_LIBRARY}.strained.tmp
else
    : > ${DATA_DIR}/${RM_LIBRARY}.strained.tmp
fi
echo "Compiling library"
if [ -f ${DATA_DIR}/trf/${RM_LIBRARY}.satellites ]; then
    cat ${DATA_DIR}/trf/${RM_LIBRARY}.satellites >> ${DATA_DIR}/${RM_LIBRARY}.strained.tmp
fi
if [ ! -s ${DATA_DIR}/${RM_LIBRARY}.strained.tmp ]; then
  echo "ERROR: TEstrainer produced an empty library" >&2
  exit 1
fi
mv ${DATA_DIR}/${RM_LIBRARY}.strained.tmp ${DATA_DIR}/${RM_LIBRARY}.strained

# Preserve intermediates until all stages have completed successfully.
echo "Removing temporary files"
rm -r ${DATA_DIR}/*/split/
rm -r ${DATA_DIR}/tmp/
find ${DATA_DIR}/ -mindepth 1 -name "run_*" -exec rm -r {} +
if [[ $RUNS -gt 0 ]]; then rm ${GENOME}.n*; fi
