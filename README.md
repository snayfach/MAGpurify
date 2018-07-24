# MAG-purify
Decontamination of (metagenome-assembled) genomes

## Installation

Install 3rd party programs:  	

* [NCBI BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download), [Prodigal](https://github.com/hyattpd/Prodigal), [HMMER](http://hmmer.org/download.html), and [LAST](http://last.cbrc.jp)
* Make sure these programs are located on your PATH

Install required python libraries:
`pip install --user pandas numpy sklearn biopython`

Clone the repo from github:  
`git clone https://github.com/snayfach/MAG-purify`

<b>Download the reference database:</b>  

From your browser:   
[http://lighthouse.ucsf.edu/IGGdb/mag-purify\_db_v1.0.tar.gz](http://lighthouse.ucsf.edu/IGGdb/mag-purify_db_v1.0.tar.gz)

Or, from the command line:   
on Unix: `wget http://lighthouse.ucsf.edu/IGGdb/mag-purify_db_v1.0.tar.gz`  
on OSX: `curl http://lighthouse.ucsf.edu/IGGdb/mag-purify_db_v1.0.tar.gz > mag-purify_db_v1.0.tar.gz`

Unpack the database: `tar -zxvf mag-purify_db_v1.0.tar.gz`   



### Testing

To list options, use: `python run_qc.py -h`  

Run MAG-cleaner on a very small test dataset
`python run_qc.py test.fna test -d mag-purify_db_v1.0`


