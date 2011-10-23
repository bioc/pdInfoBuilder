\name{NEWS}
\title{News for Package 'pdInfoBuilder'}

\section{Changes in version 1.17}{
  \subsection{USER VISIBLE CHANGES}{
    \itemize{
      \item The annotation packages for SNP chips no longer have a
      'fragment_length*' column in the featureSet/featureSetCNV tables.
      Fragment length now has its own table, 'fragmentLength' (for SNP probes)
      and 'fragmentLengthCNV' (for CN probes).
      \item Rda files available through annotation packages created by
      pdInfoBuilder > 1.17 are compressed using 'xz'.
      \item Fixed License field on the template.
      \item Added an additional template for SNP arrays supported by CRLMM.
    }
  }
}

\section{Changes in version 1.15}{
  \subsection{USER VISIBLE CHANGES}{
    \itemize{
      \item Added support to Axiom chips: 'AffySNPPDInfoPkgSeed2' gets a
      new slot 'axiom', which takes TRUE or FALSE, depending on the
      specific array type.
    }
  }
  
}

\section{Changes in version 1.13}{
  \subsection{BUG FIXES}{
    \itemize{
      \item Added 'ZipData: no' to the template, which
        caused problems on Windows systems.
    }
  }
}

\section{Changes in version 1.11}{
  
  \subsection{USER VISIBLE CHANGES}{
    \itemize{
      \item For Gene and Exon arrays, all PM probes are in
        the pmfeature table (matching what happens with
        expression arrays)
      \item Updated package template to require oligo (>= 1.11.18)
    }
  }

  \subsection{BUG FIXES}{
    \itemize{
      \item Fixed creation of indices for Gene arrays
    }
  }
}
