\NeedsTeXFormat{LaTeX2e}[1999/12/01]
\ProvidesClass{poster}
[2026/11/06 v1.2.0 academic poster class]
% options declare
\RequirePackage{kvoptions}
\SetupKeyvalOptions{family=ap,prefix=ap@}
\DeclareStringOption[ffffff]{cbg}
\DeclareStringOption[000000]{cfg}
\DeclareStringOption[0099cc]{cln}
\DeclareStringOption[1189mm]{width}
\DeclareStringOption[841mm]{height}
\DeclareStringOption[12mm]{fs}
\DeclareStringOption[6mm]{gap}
\DeclareStringOption[2mm]{sep}
\DeclareStringOption[2mm]{lw}
\DeclareStringOption[3]{ncol}
\DeclareStringOption[numeric-comp]{bibstyle}
% options process
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{article}}
\ExecuteOptions{}
\ProcessKeyvalOptions*
\ProcessOptions\relax
\LoadClass{article}
% colors
\RequirePackage{xcolor}
\definecolor{cbg}{HTML}{\ap@cbg}
\definecolor{cfg}{HTML}{\ap@cfg}
\definecolor{cln}{HTML}{\ap@cln}
\pagecolor{cbg}
% fonts
\RequirePackage[fontsize=\ap@fs]{fontsize}
\newcommand{\titlefont}{\huge}
\newcommand{\secfont}{\large\bfseries}
\newcommand{\parfont}{\bfseries}
\newcommand{\itemfont}{\color{cfg}\bfseries}
% lengths
\newlength{\gap}\setlength{\gap}{\ap@gap}
\newlength{\sep}\setlength{\sep}{\ap@sep}
\setlength{\textfloatsep}{\sep}
\setlength{\intextsep}{\sep}
\setlength{\parindent}{0pt}
\setlength{\parskip}{0pt}
\frenchspacing
% misc
\pagestyle{empty}
\newcommand{\hsep}{\hspace{\sep}\ignorespaces}
\renewcommand{\paragraph}[1]{\par{\parfont{#1}}\hsep}
\renewcommand{\mp}[3][t]{%
  \minipage[#1]{\dimexpr#2\textwidth}\vskip0pt#3\par
  \endminipage\ignorespaces}
\renewcommand{\wp}[3][t]{\hsep\mp[#1]{-\sep+#2}{#3}}
% columns
\RequirePackage{multicol}
\let\columnsep\gap
\let\multicolsep\gap
\RequirePackage[
  paperwidth=\ap@width,
  paperheight=\ap@height,
  margin=\gap,
]{geometry}
% tcolorbox
\RequirePackage{tcolorbox}
\tcbuselibrary{breakable}
\tcbset{
  top=\sep,
  left=\sep,
  right=\sep,
  bottom=\sep,
  before skip=\gap,
  after skip=\gap,
  arc=.5\gap,
  boxrule=\ap@lw,
  boxsep=0pt,
  parbox=false,
  fonttitle=\secfont,
  colframe=cfg,
  colback=white}
\newif\ifintcb
\newcommand{\esctcb}{\ifintcb\endtcolorbox\intcbfalse\fi}
% sections
\renewcommand{\section}[2][]{\esctcb%
  \tcolorbox[#1,adjusted title=#2]\intcbtrue}
\let\ocolumnbreak\columnbreak
\def\columnbreak{\esctcb\ocolumnbreak}
\AtEndDocument{\esctcb\end{multicols*}}
% title block
\RequirePackage{authblk}
\setlength{\affilsep}{0pt}
\renewcommand\AB@affilsepx{\hsep}
\renewcommand{\maketitle}[1][t]{
  \tcolorbox
  \mp[#1]{\w@tbl+0}{\@tbl}
  \mp[#1]{-\w@tbl-\w@tbr+1}{{\titlefont{\@title}}\vskip\sep{\@author}}
  \mp[#1]{\w@tbr+0}{\@tbr}
  \endtcolorbox
  \begin{multicols*}{\ap@ncol}}
\newlength{\w@tbl}\newlength{\w@tbr}
\newcommand{\tbl}[1]{\settowidth\w@tbl{#1}\def\@tbl{#1}}\tbl{}
\newcommand{\tbr}[1]{\settowidth\w@tbr{#1}\def\@tbr{#1}}\tbr{}
% floats
\RequirePackage{float,subcaption}
\captionsetup{
  skip=\sep,
  font={small},
  labelfont={bf,color=cfg},
  singlelinecheck=false}
\renewcommand{\fps@figure}{H}
\renewcommand{\fps@table}{H}
% enumitem
\RequirePackage{enumitem}
\setlist{
  labelsep=\sep,
  itemsep=0pt,
  topsep=0pt,
  parsep=0pt,
  leftmargin=*}
\setlist[itemize]{label=\itemfont{\textbullet}}
\setlist[enumerate]{label=\itemfont{\arabic*.}}
% bib
\RequirePackage[
  natbib=true,
  style=\ap@bibstyle,
  sorting=none,
  isbn=false,
  date=year,
]{biblatex}
\let\cite\supercite
\AtEveryBibitem{\iffieldundef{doi}{}{\clearfield{url}}}
\DeclareFieldFormat{labelnumberwidth}{{\itemfont{#1\adddot}}}
\renewcommand{\bibfont}{\scriptsize}
\let\oprintbibliography\printbibliography
\renewcommand\printbibliography{%
  \section{References}\oprintbibliography[heading=none]}
% hyperref
\RequirePackage[
  colorlinks,
  linkcolor=cln,
  citecolor=cln,
  urlcolor=cln,
]{hyperref}
