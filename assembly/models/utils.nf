def get_software(String software_dir, String conda_path) {
    python_env = "${conda_path}/python"
    return [
        samtools: "${software_dir}/bin/samtools",
        guppy_barcoder: "${software_dir}/bin/guppy_barcoder",
        bcftools: "${software_dir}/bin/bcftools",
        filtlong: "${software_dir}/Filtlong/bin/filtlong",
        kraken2: "${software_dir}/kraken-2.12/kraken2",
        minimap2: "${software_dir}/minimap2/minimap2",
        python: "${python_env}/bin/python",
        seqtk : "/usr/bin/seqtk",
        blastn: "${software_dir}/ncbi-blast-2.14.0+/bin/blastn"
    ]
}

def q2p(int qual){
    return (1-10**(-qual/10)) * 100
}

def print_error(String message) {
    String line = "=" * message.length()
    println("\033[31m${line}\n${message}\n${line}\033[0m")
    System.exit(1)
}