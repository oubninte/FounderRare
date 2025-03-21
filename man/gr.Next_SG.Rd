% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gr.Next_SG.R
\name{gr.Next_SG}
\alias{gr.Next_SG}
\title{Construct the Next SG range for a particular SG.}
\usage{
gr.Next_SG(gr.SG, minL = 400, Overlap = 0)
}
\arguments{
\item{gr.SG}{Grange of concerned SG}

\item{minL}{Minimum length of SG Grange in Kbp, 400Kbp by default}

\item{Overlap}{Overlap between current SG and next SG in bp; equal to 0 by default}
}
\value{
Grange of next SG
}
\description{
Construct the Next SG range for a particular SG.
}
\examples{
gr.Next_SG(gr.SG)
gr.Next_SG(gr.SG, minL=400)
gr.Next_SG(gr.SG, minL=400, Overlap=0)


}
