\name{NEWS}
\title{News for Package 'pdInfoBuilder'}

\section{Changes in version 1.32}{
  \subsection{BUG FIXES}{
    \itemize{
      \item Fixed dependency changes for RSQLite 1.0.0
     }
   }
   \subsection{USER VISIBLE CHANGES}{
      \itemize{
         \item GenericArray is now available to accommodate oligo support for custom CDFs
      }
    }
}

\section{Changes in version 1.28}{
  \subsection{USER VISIBLE CHANGES}{
    \itemize{
      \item Changed dependency on RSQLite. Now, version >= 1.0.0 is required.
      \item GenericArray is now available to support customized annotations.
    }
  }
}

\section{Changes in version 1.22}{
  \subsection{BUG FIXES}{
    \itemize{
      \item Removed dependency on RConverters.h (Prof Ripley)
    }
  }
}

\section{Changes in version 1.20}{
  \subsection{USER VISIBLE CHANGES}{
    \itemize{
      \item The annotation packages for Affymetrix Tiling arrays now
      contain strand information (suggestion and implementation: Kristof De Beuf)
    }
  }
}

\section{Changes in version 1.18}{
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
