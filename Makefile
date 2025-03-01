# MakeFile
install:
	Rscript -e "source('install.R')"
	Rscript -e "if (!require('devtools' quietly = TRUE)){
     	install.packages('devtools')}
	devtools::install_github('mamouzgar/hafez')"


