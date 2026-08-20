\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{ms}
[2024/02/26 v0.1.4 Manuscript Class]
% options
\RequirePackage{kvoptions}
\SetupKeyvalOptions{family=ms,prefix=ms@}
\DeclareStringOption[numeric-comp]{bibstyle}
\DeclareStringOption[35mm]{margin}
\DeclareStringOption[onehalf]{spacing}
\DeclareBoolOption[false]{blind}
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{article}}
\ExecuteOptions{}
\ProcessKeyvalOptions*
\ProcessOptions\relax
\LoadClass{article}
% spacing
\RequirePackage[margin=\ms@margin]{geometry}
\setlength{\parskip}{1ex}
\setlength{\parindent}{0ex}
\setlength{\abovecaptionskip}{\medskipamount}
\setlength{\belowcaptionskip}{\medskipamount}
\newlength{\metalabel}\setlength{\metalabel}{.15\linewidth}
\RequirePackage[\ms@spacing spacing]{setspace}
% meta-data
\RequirePackage{authblk}
\newcommand{\blind}[1]{\ifms@blind\textsc{redacted}\else#1\fi}
\def\meta@block{}
\newcommand{\meta@part}[2]{%
  \hspace{\metalabel}\llap{\textbf{#1}}\hspace{1em}%
  \parbox[t]{\dimexpr\linewidth-\metalabel-1em}{#2}\par\medskip}
\newcommand{\metafield}{\@ifstar\blind@metafield\open@metafield}
\newcommand{\open@metafield}[3][]{%
  \g@addto@macro{\meta@block}{%
    \meta@part{\ifempty{#1}{\MakeUppercase#2}{#1}}{#3}}
  \expandafter\gdef\csname @#2\endcsname{#3}}
\newcommand{\blind@metafield}[3][]{\open@metafield[#1]{#2}{\blind{#3}}}
\renewcommand{\title}[2][]{%
  \metafield{title}{#2}%
  \ifempty{#1}{}{\metafield[Short Title]{shorttitle}{#1}}
  \g@addto@macro{\meta@block}{\meta@part{Authors}{\blind{\@author}}\medskip}}
\renewcommand{\date}[1]{\metafield{date}{#1}}
\renewcommand{\abstract}[1]{\metafield{abstract}{#1}}
\newcommand{\metabreak}{\g@addto@macro{\meta@block}{\clearpage}}
% maketitle
\renewcommand{\maketitle}{%
  \clearpage\setcounter{page}{1}%
  \if@twocolumn\onecolumn\meta@block\twocolumn
  \else\meta@block\fi
  \clearpage}
% appendix - TODO: option for part vs section
\let\xappendix\appendix
\renewcommand{\appendix}{%
  \clearpage\setcounter{page}{1}%
  \meta@part{}{\textbf{\appendixname}}
  \meta@part{Title}{\@title}
  \meta@part{Authors}{\blind{\@author}}\medskip
  \meta@part{Date}{\@date}
  \xappendix
  \renewcommand{\partname}{Appendix}
  \renewcommand{\thepart}{\Alph{part}}
  \renewcommand{\thesection}{\thepart.\arabic{section}}
  \refsection
  \counterwithin{section}{part}
  \counterwithin{table}{part}
  \counterwithin{figure}{part}
  \counterwithin{equation}{part}
  \clearpage}
% cover
\newcommand{\address}[3][Submission of a manuscript to \emph{\@journal}]{%
  \thispagestyle{empty}%
  \@date\medskip\par
  \parbox[t]{.5\linewidth}{\raggedright#2}%
  \parbox[t]{.5\linewidth}{\raggedleft#3}\bigskip\par
  \textbf{Re: #1}\medskip\par}
% bib
\RequirePackage[
  natbib=true,
  style=\ms@bibstyle,
  sorting=none,
  url=true,
  doi=false,
  isbn=false,
  date=year,
]{biblatex}
% helpers
\renewcommand{\paragraph}[1]{\par\textbf{#1}\hskip.5em}
\newcommand{\ifempty}[3]{\if\relax\detokenize{#1}\relax#2\else#3\fi}
