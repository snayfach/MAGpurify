# MAG-purify
Decontamination of (metagenome-assembled) genomes

## Installation

Install 3rd party programs:  	

* [NCBI BLAST] (ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/), [Prodigal] (https://github.com/hyattpd/Prodigal), [HMMER] (http://hmmer.org/download.html), and [LAST] (http://last.cbrc.jp)
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
`python run_qc.py ... --verbose`



# MAG cleaner

## Dependencies

1. Python 2
2. Python modules: pandas, numpy, sklearn, biopython
3. External programs:
	

## Database

The iMAGen database can be downloaded from dropbox via [http] (https://www.dropbox.com/sh/xye5c3lpm1oscfp/AAC8Wf6UOMnUQDNLVqCKvKDea?dl=0),

Or using wget:
`wget `

After downloading, add the directory to your environment:  
`export IMAGEN_DB=/path/to/iMAGen-db`
