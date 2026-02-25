import os
import sys
import subprocess
import urllib.request
from pathlib import Path
import yaml

def running_tea(stage="Starting Earl Grey"):
    """Display ASCII art tea cup with stage name"""
    tea_art = rf"""    
          )  (
         (   ) )
         ) ( (
       _______)_
    .-'---------|  
   ( C|/\/\/\/\/|
    '-./\/\/\/\/|
      '_________'
       '-------'
    <<< {stage} >>>"""
    print(tea_art)

def convert_seconds(seconds):
    """Convert seconds to HH:MM:SS.ss format"""
    seconds = int(seconds)
    h = seconds // 3600
    m = seconds % 3600 // 60
    s = seconds % 60
    return f"{h:02d}:{m:02d}:{s:05.2f}"


# --------------------------------------------------
# Message helpers (consistent style)
# --------------------------------------------------
def msg_header(title):
    print("\n" + "=" * 60)
    print(f"{title}")
    print("=" * 60)

def msg_info(text):
    print(f"[INFO] {text}")

def msg_warn(text):
    print(f"[WARNING] {text}")

def msg_error(text):
    sys.exit(f"\n[ERROR] {text}\n")

# --------------------------------------------------
# Validation + defaults with custom messages
# --------------------------------------------------
def validate_parameters(config, defaults=None, outfile = None):

    msg_header("Earl Grey configuration check")

    # ---- Required parameters ----
    required = ['genome', 'species', 'output_dir']
    for param in required:
        if not config.get(param):
            msg_error(f"Required parameter '{param}' not specified in config file")

    msg_info("Required parameters detected")

    # ---- Defaults + messages ----
    defaults = {
        'threads': (1, "{} cores will be used"),
        'num': (10, "De Novo Sequences will be extended through a maximum of {} iterations"),
        'no_seq': (20, "{} sequences will be used in BEAT consensus generation"),
        'cluster': ('no', None),
        'softMask': ('no', None),
        'margin': ('no', None),
        'Flank': (1000, "Blast, extend, align, trim process will add {}bp to each end in each iteration"),
        'min_seq': (3, "Blast, extend, align, trim process will require {} sequences to generate a new consensus sequence"),
        'run_heliano': ('FALSE', None)
    }

    msg_header("Parameter values")

    for param, (default_val, message_template) in defaults.items():
        if not config.get(param):
            config[param] = default_val
            if message_template:
                msg_info(message_template.format(default_val))
        else:
            if message_template:
                msg_info(message_template.format(config[param]))

    # --------------------------------------------------
    # Verbose user messages 
    # --------------------------------------------------
    msg_header("Pipeline behaviour")

    # RepeatMasker
    if not config.get('RepSpec') and not config.get('startCust'):
        msg_info("RepeatMasker species not specified, running Earl Grey without an initial mask with known repeats")
    else:
        msg_info("Running with an initial mask using known repeats")

    # Clustering
    if config['cluster'] == 'yes':
        msg_warn("TE consensus sequences will be clustered (may affect subfamilies and create chimeras)")
    else:
        msg_info("TE consensus sequences will not be clustered")

    # SoftMask
    if config['softMask'] == 'yes':
        msg_info("Softmasked genome will be generated")
    else:
        msg_info("Softmasked genome will not be generated")

    # Margin
    if config['margin'] == 'yes':
        msg_info("Short TE sequences (<100bp) will be removed")
    else:
        msg_info("Short TE sequences (<100bp) will not be removed")

    # Helitrons
    if config['run_heliano'] == 'yes':
        msg_info("HELITRON detection will be run using HELIANO")
    else:
        msg_info("HELITRON detection will not be run")

    # --------------------------------------------------
    # Structural validation 
    # --------------------------------------------------
    msg_header("Input validation")

    if not isinstance(config['genome'], dict):
        msg_error("'genome' must be a dictionary: species → fasta")

    for sp, path in config['genome'].items():
        if not Path(path).exists():
            msg_error(f"Genome file for species '{sp}' not found: {path}")

    if not isinstance(config['species'], list):
        msg_error("'species' must be a list")

    for sp in config['species']:
        if sp not in config['genome']:
            msg_error(f"Species '{sp}' listed but no genome provided")

    if not os.path.isdir(config["script_dir"]):
        msg_error("ERROR: Script directory variable not set, please run the configure script in the Earl Grey directory before attempting to run Earl Grey")
    
    testrainer_dir = os.path.join(config["script_dir"], "TEstrainer")
    if not os.path.isdir(testrainer_dir):
        msg_error("ERROR: teStrainer module not found, please check all modules are present and run the configure script in the Earl Grey directory before attempting to run Earl Grey")

    # --------------------------------------------------
    # Output setup
    # --------------------------------------------------
    msg_header("Output setup")

    outdir = Path(config['output_dir'])
    outdir.mkdir(parents=True, exist_ok=True)
    msg_info(f"Output directory: {outdir}")

    for sp in config['species']:
        sp_dir = outdir / f"{sp}_EarlGrey"
        sp_dir.mkdir(parents=True, exist_ok=True)
        msg_info(f"Created directory: {sp_dir}")

    # Save validated config (optional but recommended)
    if outfile:
        with open(outfile, "w") as f:
            yaml.safe_dump(config, f, sort_keys=False)
            msg_info(f"Validated configuration saved to: {outfile}")


    msg_header("Configuration complete")
    msg_info("All checks passed. Workflow ready to start.\n")
    
    print("\nPlease cite the following paper when using this software:")
    print("Baril, T., Galbraith, J. and Hayward, A., 2024. Earl Grey: a fully automated user-friendly transposable element annotation and analysis pipeline. Molecular Biology and Evolution, 41(4), p.msae068. \n")
    

    return config


def make_directories(directory, species, RepSpec=None, startCust=None, run_heliano=None):
    outdir = os.path.join(directory, f"{species}_EarlGrey")
    os.makedirs(outdir, exist_ok=True)

    if RepSpec or startCust:
        os.makedirs(os.path.join(outdir, f"{species}_RepeatMasker"), exist_ok=True)
    os.makedirs(os.path.join(outdir, f"{species}_Database"), exist_ok=True)
    os.makedirs(os.path.join(outdir, f"{species}_RepeatModeler"), exist_ok=True)
    os.makedirs(os.path.join(outdir, f"{species}_strainer"), exist_ok=True)
    os.makedirs(os.path.join(outdir, f"{species}_RepeatMasker_Against_Custom_Library"), exist_ok=True)
    os.makedirs(os.path.join(outdir, f"{species}_RepeatLandscape"), exist_ok=True)
    os.makedirs(os.path.join(outdir, f"{species}_mergedRepeats"), exist_ok=True)
    os.makedirs(os.path.join(outdir, f"{species}_summaryFiles"), exist_ok=True)
    if run_heliano:
        os.makedirs(os.path.join(outdir, f"{species}_heliano"), exist_ok=True)
    return outdir

#CHECK WITH TOBY: Should this be in the pangenome pipeline? 
# Or is this part of what gets set up by first running EarlGrey normally?
def check_biocontainer():
    """Check and configure biocontainer installation if needed"""
    repeatmasker_lib = "/usr/local/share/RepeatMasker/Libraries/"
    dfam_file = os.path.join(repeatmasker_lib, "Dfam.h5")
    
    if os.path.isdir(repeatmasker_lib) and os.path.isfile(dfam_file):
        try:
            with open(dfam_file, 'r') as f:
                content = f.read()
                if "Placeholder" in content:
                    response = input("Are you using the biocontainer installation? If yes, we can configure RepeatMasker to work for you:[YyNn]")
                    if response.lower() in ['y', 'yes']:
                        configure_biocontainer(repeatmasker_lib)
                    else:
                        sys.exit(0)
        except Exception as e:
            print(f"Warning: Could not read Dfam.h5 file: {e}")

def configure_biocontainer(lib_path):
    """Configure biocontainer RepeatMasker installation"""
    try:
        os.chdir(lib_path)
        
        # Download Dfam database
        print("Downloading Dfam database...")
        urllib.request.urlretrieve(
            "https://dfam.org/releases/Dfam_3.7/families/Dfam_curatedonly.h5.gz",
            "Dfam_curatedonly.h5.gz"
        )
        
        # Extract and replace
        subprocess.run(["gunzip", "Dfam_curatedonly.h5.gz"], check=True)
        if os.path.exists("Dfam.h5"):
            os.rename("Dfam.h5", "Dfam.h5.bak")
        os.rename("Dfam_curatedonly.h5", "Dfam.h5")
        
        # Configure RepeatMasker
        os.chdir("/usr/local/share/RepeatMasker/")
        subprocess.run([
            "perl", "configure",
            "-rmblast_dir=/usr/local/bin",
            "-libdir=/usr/local/share/RepeatMasker/Libraries",
            "-trf_prgm=/usr/local/bin/trf",
            "-default_search_engine=rmblast"
        ], check=True)
        
        print("Biocontainer configuration completed successfully")
        
    except Exception as e:
        print(f"Error configuring biocontainer: {e}")
        sys.exit(1)


#CHECK WITH TOBY: Should this be in the pangenome pipeline? 
# Or is this part of what gets set up by first running EarlGrey normally?
def check_dfam39():
    """Check for Dfam 3.9 configuration requirements"""
    try:
        # Get RepeatMasker path
        rm_path = subprocess.run(["which", "RepeatMasker"], 
                                capture_output=True, text=True, check=True).stdout.strip()
        library_path = rm_path.replace("/bin/RepeatMasker", "/share/RepeatMasker/Libraries/famdb")
        
        expected_file = os.path.join(library_path, "dfam39_full.0.h5")
        config_file = os.path.join(library_path, "rmlib.config")
        complete_marker = os.path.join(library_path, ".earlgrey.config.complete")
        
        # Check if configuration is needed
        #Check with Toby should these be "OR" instead of "AND"? Should there be a more complex check?
        if (os.path.isdir(library_path) and 
            len([f for f in os.listdir(library_path) if os.path.isfile(os.path.join(library_path, f))]) == 2 and
            os.path.isfile(expected_file) and 
            os.path.isfile(config_file) and
            not os.path.exists(complete_marker)):
            
            generate_dfam39_config_script(library_path, rm_path)
            sys.exit(2)
        #Check with Toby for the following
        #else :
        #    open(complete_marker, "w").close()
        return complete_marker
            
    except subprocess.CalledProcessError:
        print("Warning: Could not locate RepeatMasker installation")
    except Exception as e:
        print(f"Warning: Error checking Dfam 3.9 configuration: {e}")

#CHECK WITH TOBY: Should this be in the pangenome pipeline? 
# Or is this part of what gets set up by first running EarlGrey normally?
def generate_dfam39_config_script(library_path, rm_path):
    """Generate configuration script for Dfam 3.9"""
    script_path = os.path.join(os.getcwd(), "configure_dfam39.sh")
    
    print("WARNING: Earl Grey v6.1.0 has updated to Dfam v3.9.")
    print("Before using Earl Grey, you MUST download the required partitions from Dfam")
    print(f"Configuration script generated: {script_path}")
    
    script_content = f"""#!/bin/bash
    # Configuration script for Dfam 3.9

    # Change directory to the famdb library location
    cd {library_path}/

    # Download partitions (modify the range [0-16] as needed)
    curl -o 'dfam39_full.#1.h5.gz' 'https://dfam.org/releases/current/families/FamDB/dfam39_full.[0-16].h5.gz'

    # Decompress files
    gunzip *.gz

    # Move to RepeatMasker directory
    cd {rm_path.replace('/bin/RepeatMasker', '/share/RepeatMasker/')}

    # Backup existing files
    mv {library_path}/min_init.0.h5 {library_path}/min_init.0.h5.bak

    # Reconfigure RepeatMasker
    perl ./configure \\
        -libdir {library_path.replace('/famdb', '')} \\
        -trf_prgm {rm_path.replace('/bin/RepeatMasker', '/bin/trf')} \\
        -rmblast_dir {rm_path.replace('/bin/RepeatMasker', '/bin')} \\
        -hmmer_dir {rm_path.replace('/bin/RepeatMasker', '/bin')} \\
        -abblast_dir {rm_path.replace('/bin/RepeatMasker', '/bin')} \\
        -crossmatch_dir {rm_path.replace('/bin/RepeatMasker', '/bin')} \\
        -default_search_engine rmblast

    # Mark configuration as complete
    touch {library_path}/.earlgrey.config.complete
    """
    
    with open(script_path, 'w') as f:
        f.write(script_content)
    
    #os.chmod(script_path, 0o755)
    print(f"Make the script executable and run: chmod +x {script_path} && ./{script_path}")

