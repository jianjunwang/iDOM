% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DOM.H2.R
\name{DOM.H2}
\alias{DOM.H2}
\title{The specialization index H2' of DOM-microbe associations}
\usage{
DOM.H2(
  Comm.Microbe,
  Comm.DOM,
  occurrence.threshold = 0.5,
  threshold.r = 0.3,
  N = 100,
  Null.model = "swap.web"
)
}
\arguments{
\item{Comm.Microbe}{Sample/Microbe matrix with samples in the rows and Microbial species in the columns (compositional/abundant data).}

\item{Comm.DOM}{Sample/DOM matrix with samples in the rows and DOM molecules in the columns (compositional/abundant data).}

\item{occurrence.threshold}{The threshold of retaining bacterial species or DOM molecules observed in more than half of the total samples. Default: 0.5.}

\item{threshold.r}{The threshold of including correlations between DOM molecules and bacterial species in bipartite networks. Defaulted: 0.3.}

\item{N}{Number of null models to be generated; defaults to 100 (more might be better, less probably not).}

\item{Null.model}{Null model type. Can be given as an integer or name: 1/"r2dtable", 2/"swap.web", 3/"vaznull", 4/"shuffle.web"; allows for partial match of names.}
}
\value{
Returns a data.frame, which contains standardised H2', Observed H2', P-value and Network type. Standardised: Standardised H2'. It is standardised by using a null modelling approach.
Observed: Observed H2'. It ranges between 0 (complete generalization) and 1 (complete specialization).p.value: The significance of differences between observed and random H2'.Network.type: Network type. Can be given as a name: Full (full network),
Negative (negative network), Positive (positive network).
}
\description{
This function calculates the network-level specialization of all interacting trophic levels in DOM-microbe bipartite networks, including full, negative and positive networks.
}
\details{
H2' is a network-level property to describe how much the two trophic
levels are interacting with each other in a bipartite network. For example,
H2' is used to quantify the specialization of DOM-microbe associations at a
network level. Specifically, elevated H2' values convey that there is a high
degree of specialization between DOM and microbes. By contrast, lower H2'
values reflect a more generalized bipartite network where different DOM
molecules can be used by a large range of bacterial taxa.
}
\examples{
\dontrun{
# Example data of a Microbial compositional table (50 samples by 100 bacterial species)
Microbial.data
# Example data of a DOM compositional table (50 samples by 100 DOM molecules)
DOM.data
# Calculation of H2' index
DOM.H2(Comm.Microbe = Microbial.data,
       Comm.DOM = DOM.data,
       occurrence.threshold = 0.5,
       threshold.r = 0.3,
       N = 100,
       Null.model = "swap.web")
}
}
\references{
Hu A, Choi M, Tanentzap AJ, Liu J, Jang K-S, Lennon JT, Liu Y, Soininen J, Lu X, Zhang Y, Shen J, Wang J.
Quantifying microbial associations of dissolved organic matter under global change. \emph{bioRxiv}, 2021.

Hu A, Choi M, Tanentzap AJ, Liu J, Jang K-S, Lennon JT, Liu Y, Soininen J, Lu X, Zhang Y, Shen J, Wang J.
Quantifying the associations between dissolved organic matter and microbes
under global change. \emph{Nature communications}, 2022.
}
\seealso{
\code{\link[SpiecEasi]{sparcc}}
\code{\link[bipartite]{nullmodel}}, \code{\link[bipartite]{networklevel}}
}
