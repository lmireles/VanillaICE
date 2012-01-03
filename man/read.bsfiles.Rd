\name{read.bsfiles}
\alias{read.bsfiles}
\title{
  Read BeadStudio/GenomeStudio processed data.
}
\description{
  Read BeadStudio/GenomeStudio processed data and return an array of log
  R ratios and B allele frequencies.
}
\usage{
read.bsfiles(path = "", filenames, ext = "", row.names = 1, sep = "\t", lrr.colname = "Log.R.Ratio", baf.colname = "B.Allele", drop = FALSE, colClasses, nrows = 1800000, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{path}{
    character: path to plain text files containing BeadStudio processed
    data
}
  \item{filenames}{
    character: name of file(s)
}
  \item{ext}{
    character: filename extension
}
\item{row.names}{
  As in read.table. By default, the first column is assumed to be the
  feature identifiers.
}
  \item{sep}{
    As in read.table.
  }

  \item{lrr.colname}{

    character: used to grep for the log R ratios in the header. E.g.,
    \code{grep(lrr.colname, header) } should return a length 1 vector,
    where header is a vector of the column labels.

  }

  \item{baf.colname}{

    character: used to grep for the B allele frequency in the
    header. E.g., \code{grep(baf.colname, header) } should return a
    length 1 vector, where header is a vector of the column labels.

  }

  \item{drop}{
      Logical:  if TRUE, dimnames will not be returned
  }

  \item{colClasses}{

    Vector as in \code{read.table}.  Note that if colClasses is not
    specified, the colClasses will be defined by reading in the first
    few rows.  "NULL" will be assigned to all columns not containing B
    allele frequencies or log R ratios.

  }

  \item{nrows}{

    As in \code{read.table}.

  }

  \item{\dots}{

    Additional arguments passed to \code{read.table}.

}
}

%\details{
%
%
%}

\value{

  A 3 dimensional array: features x statistic (lrr or baf) x sample

}

\author{
R. Scharpf
}


\seealso{
  \code{\link{read.table}}
}
\examples{
path <- system.file("extdata", package="VanillaICE")
filename <- list.files(path, pattern="LRRandBAF", full.names=TRUE)
dat <- read.bsfiles(filenames=filename)
}
\keyword{IO}

