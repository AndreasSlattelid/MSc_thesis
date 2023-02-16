%% Identification part:

\NeedsTeXFormat{LaTeX2e}[1995/12/01]  % Any version, really!
\ProvidesClass{uiomasterbeta}[2022/11/29 v1.1 UiO document class for
  master's thesis by Dag Langmyhr, Martin Helsø and Karoline Moe]

%% Required packages:
\RequirePackage[T1]{fontenc}
\RequirePackage{xifthen}

%% Option declaration part:
\DeclareOption{bachelor}{\def \uiom@base {article}}
\DeclareOption{binding}{\def \uiom@bind {24pt}}
\DeclareOption{font=cmr}{\def \uiom@font {\RequirePackage{lmodern}
    \RequirePackage[scaled]{helvet}}}
\DeclareOption{font=garamond}{\def \uiom@font {\RequirePackage[garamond]{mathdesign}
    \RequirePackage[scaled]{helvet}}}
\DeclareOption{font=times}{\RequirePackage{txfonts}}
\DeclareOption{master}{\def \uiom@base {report}}

\def \uiom@bind {0pt}

\DeclareOption*{\PassOptionsToClass{\CurrentOption}{\uiom@base}}

\ExecuteOptions{font=cmr, master}
\ProcessOptions \relax

%% Base class load part:
\ifthenelse{\equal{\uiom@base}{report}}
  {\LoadClass[a4paper,11pt,twoside,openright]{report}}
  {\LoadClass[a4paper,11pt]{article}}

%% Margins:
\RequirePackage[vmargin=3cm,bindingoffset=\uiom@bind]{geometry}

%% Text adjustments:
\pretolerance = 3000
\tolerance = 6000   \hbadness = \tolerance
\setlength{\parskip}{0cm plus 1mm}

%% Select a font
\uiom@font

%% Chapter and section headers:
\RequirePackage[]{titlesec}
\titleformat{\part}[display]
  {\sf \bfseries \Large \centering}{Part \thepart}{15pt}{\Huge\sf}
\titleformat{\chapter}[display]
  {\sf\LARGE}{\chaptertitlename\ \thechapter}{15pt}{\Huge \sf}
\titleformat{\section}
  {\sf\Large}{\thesection}{1em}{}
\titleformat{\subsection}
  {\sf\large}{\thesubsection}{1em}{}
\titleformat{\subsubsection}
  {\sf\normalsize}{\thesubsubsection}{1em}{}
\titleformat{\paragraph}[runin]
  {\sf\normalsize}{\theparagraph}{1em}{}
\titleformat{\subparagraph}[runin]
  {\sf\normalsize}{\thesubparagraph}{1em}{}
\titlespacing*{\chapter} {0pt}{50pt}{40pt}
\titlespacing*{\section} {0pt}{3.5ex plus 1ex minus .2ex}{2.3ex plus .2ex}
\titlespacing*{\subsection} {0pt}{3.25ex plus 1ex minus .2ex}{1.5ex plus .2ex}
\titlespacing*{\subsubsection}{0pt}{3.25ex plus 1ex minus .2ex}{1.5ex plus .2ex}
\titlespacing*{\paragraph} {0pt}{3.25ex plus 1ex minus .2ex}{1em}
\titlespacing*{\subparagraph} {\parindent}{3.25ex plus 1ex minus .2ex}{1em}

%% Table of contents, list of figures, list of tables
\usepackage{tocloft,lipsum,pgffor}
\setcounter{tocdepth}{1}% Include up to \section in ToC
\renewcommand{\cfttoctitlefont}{\sf\LARGE}
\renewcommand{\cftpartfont}{\normalfont\sffamily \bfseries}% \part font in ToC
\renewcommand{\cftpartpagefont}{\normalfont\sffamily \bfseries}% \part font in ToC
\renewcommand{\cftchapfont}{\normalfont\sffamily}    % \chapter font in ToC
\renewcommand{\cftchappagefont}{\normalfont\sffamily}    % \chapter font in ToC
\renewcommand{\cftsecfont}{\normalfont\sf} 
\renewcommand{\cftsecpagefont}{\normalfont\sf} 
\renewcommand\cftloftitlefont{\sf\LARGE}
\renewcommand\cftlottitlefont{\sf\LARGE}
\renewcommand{\cftfigdotsep}{\cftnodots} % no dots in list of figures
\renewcommand{\cfttabdotsep}{\cftnodots} % no dots in list of tables
\renewcommand{\cftdotsep}{\cftnodots} % no dots in list of tables
\renewcommand\cftfigfont{\normalfont\sf}
\renewcommand\cfttabfont{\normalfont\sf}
\renewcommand\cftfigpagefont{\normalfont\sf}
\renewcommand\cfttabpagefont{\normalfont\sf}


%% If you want to reset chapter numbering after each new part. Needs more careful thought.
%\makeatletter
%\@addtoreset{chapter}{part}
%\makeatletter

%% If you want to use \part* for introduction and conclusion, you need to add to get TOC and header immediately after.
%\addcontentsline{toc}{part}{Introduction}
%\markboth{Introduction}{Introduction} 


%% Redefine \begin{abstract}...\end{abstract}
\renewenvironment{abstract}{%
  \ifthenelse{\equal{\uiom@base}{report}}
    {\cleardoublepage
     \markboth{\abstractname}{}%           
     \vspace*{5\baselineskip}}%
    {\bigskip}%
  \normalsize
  \begin{center}%
    {\sf \abstractname\vspace{-.5em}\vspace{\z@}}%
  \end{center}%
  \quotation}%
  \cleardoublepage
{\endquotation \bigskip}

%% Redefine \begin{sammendrag}...\end{sammendrag}
\newenvironment{sammendrag}{%
  \ifthenelse{\equal{\uiom@base}{report}}
    {\cleardoublepage
     \markboth{Sammendrag}{}%           
     \vspace*{5\baselineskip}}%
    {\bigskip}%
  \normalsize
  \begin{center}%
    {\sf Sammendrag\vspace{-.5em}\vspace{\z@}}%
  \end{center}%
  \quotation}%
  \cleardoublepage
{\endquotation \bigskip}

%% Define \begin{preface}
\newenvironment{preface}{%
  \ifthenelse{\equal{\uiom@base}{report}}
    {\chapter*{\prefacename}
      \markboth{\prefacename}{}}
    {\section*{\prefacename}}}
  {}
  
%% Define \begin{acknowledgements}
\newcommand{\acknowledgementsname}{Acknowledgements}
\newenvironment{acknowledgements}{%
  \ifthenelse{\equal{\uiom@base}{report}}
    {\chapter*{\acknowledgementsname}
      \markboth{\acknowledgementsname}{}}
    {\section*{\acknowledgementsname}}}
  {}

%% Headers and footers:
\ifthenelse{\equal{\uiom@base}{report}}
  {\RequirePackage{fancyhdr}
   \fancyhead{}
   \setlength{\headheight}{15pt}
   \fancyhead[LE]{\textsf{\small \nouppercase{\leftmark}}}
   \fancyhead[RO]{\textsf{\small \nouppercase{\rightmark}}}
   \renewcommand{\headrulewidth}{0pt}
   \fancyfoot{}
   \fancyfoot[LE,RO]{\textsf{\thepage}}
   \renewcommand{\footrulewidth}{0pt}
   \pagestyle{fancy}
   \fancypagestyle{plain}{\fancyhead{}}}
 {}

%% Figures og tables
\RequirePackage[labelfont=sf,font=small]{caption}

%% Document parts:
\newcommand{\frontmatter}{\cleardoublepage \pagenumbering{roman}}
\newcommand{\mainmatter}{\cleardoublepage \pagenumbering{arabic}}
\newcommand{\backmatter}{\cleardoublepage}

%% Lists:
\usepackage{enumitem}
\setlist[description]{font=\normalfont \sf \space}

\providecommand{\preface}{\chapter*{}\markboth{}{}}
\providecommand{\endpreface}{}

\endinput


\graphicspath{{figures/}}