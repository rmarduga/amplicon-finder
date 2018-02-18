# AmpliconFinder

Amplicon finder is a HMM-based method for detecting double minute chromosomes amplicons.

# Dependencies

AmpliconFinder requires following python libraries:

* PySam
* pomegranate
* NumPy




# Usage

<pre>
amplicon_finder.py [-h] [-t THR_NUM] [-k COV_THSH] [-i WINDOW]
                   [-g MAX_GAP] [-l MIN_AMPLICON_LENGTH]
                   [-c PATH_TO_CACHE_DIR]
                   &lttumor_bam&gt [&ltgermline_bam&gt]

positional arguments:
  &lttumor_bam&gt           Path to the tumor bam file
  &ltgermline_bam&gt        Path to the germline bam file. If provided, relative
                        coveage between tumor and germline samples will be
                        used. (Optional)

optional arguments:
  -h, --help            show this help message and exit
  -t THR_NUM, --threads THR_NUM
                        number of threads (default is 1)
  -k COV_THSH, --coverage-threshold COV_THSH
                        Coverage threshold in standard diviation units.
                        (default is 3)
  -i WINDOW, --interval WINDOW
                        Length of the interval of DNA that is being averaged.
                        It could be seen as the resolution of the method.
                        (default is 100)
  -g MAX_GAP, --max-gap MAX_GAP
                        Maximum gap in basepairs between amplicons that are
                        considered as a single amplicon. (default is 1000)
  -l MIN_AMPLICON_LENGTH, --min-amplicon-length MIN_AMPLICON_LENGTH
                        Minimum amplicon length in basepairs to be captured.
                        (default is 1000)
  -c PATH_TO_CACHE_DIR, --cache-dir PATH_TO_CACHE_DIR
                        Change path to cache dir (default is
                        './.amplicon_finder')
</pre>


# Output format
Founded amplicons will be list through the standard output in BED format.
