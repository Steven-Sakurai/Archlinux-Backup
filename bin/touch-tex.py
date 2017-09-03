 #!/usr/bin/env python 
 # -*- coding: utf-8 -*-

import sys
filename = sys.argv[1]
with open(filename, 'w') as f:
    f.write(r"""   %!TEX TS-program = xelatex
    %!TEX encoding = UTF-8 Unicode
\documentclass[a4paper,12pt]{article}
\usepackage[top=1.5in, bottom=1.5in, left=1.0in, right=1.0in]{geometry}
\usepackage{longtable}
\usepackage{geometry}
\usepackage{multirow}
\usepackage{setspace}
\usepackage{booktabs}
\usepackage{indentfirst}
\usepackage{fontspec,xltxtra,xunicode}
\usepackage{float}
\usepackage{amsmath}
\usepackage{amssymb}
\usepackage{amsthm}
\usepackage{amsfonts}
\usepackage{braket}
%\usepackage{ctex}
\usepackage[]{xeCJK}
\setCJKmainfont[BoldFont=STKaitiSC-Bold, ItalicFont=STHeitiSC-Light]{STSong}
\setCJKsansfont[BoldFont=STHeiti]{STXihei}
\setCJKmonofont{STFangsong}
%\defaultfontfeatures{Mapping=tex-text}
%\setromanfont{SimSun} %设置中文字体
\XeTeXlinebreaklocale “zh”
\XeTeXlinebreakskip = 0pt plus 1pt minus 0.1pt %文章内中文自动换行
%\newfontfamily{\H}{SimHei}
%\newfontfamily{\E}{Arial}  %设定新的字体快捷命令

%%%%%%%%%%%%%%%%%%%%%
%%%%%  MY MACRO  %%%%
%\newcommand{\lap}[1]{\mathcal{L}[{#1}]}
%\newcommand{\relap}[1]{\mathcal{L}^{-1}[{#1}]}

\title{\textbf{  Title  }}
\author{沈士恒 1500015941}
\date{\today}
\begin{document}
\maketitle


\begin{figure}[H]
\centering
\includegraphics[scale=.85]{*.png}
\caption{}
\end{figure}

\begin{longtable}{|l|r|r|}
\hline
\caption{}
\end{longtable}

\end{document}
""")
