{smcl}
{* 20 Feb 2002/19 Jul 2006}{...}
{hline}
help for {hi:xgmlt}
{hline}

{title:'Logit' Estimation with Spatial Standard Errors}

{title:Description}

{p 4 4 2}
{cmd:xgmlt} 'Logit' Estimation with Spatial Standard Errors.

{title:Examples}

{phang2} import excel "https://raw.githubusercontent.com/timbulwidodostp/xgmlt/main/xgmlt/xgmlt.xlsx", sheet("Sheet1") firstrow clear {p_end}
{phang2} sum dep {p_end}
{phang2} generate newdep = cond(dep<7.2,0,1) {p_end}
{phang2} gen const=1 {p_end}
{phang2} gen cutoff1=4 {p_end}
{phang2} gen cutoff2=4 {p_end}
{phang2} xgmlt C1 C2 cutoff1 cutoff2 newdep const indep1, xreg(2) coord(2) {p_end}
{phang2} matrix se4 = cov_nd {p_end}
{phang2} replace cutoff1 = 4 {p_end}
{phang2} replace cutoff2 = 6 {p_end}
{phang2} xgmlt C1 C2 cutoff1 cutoff2 newdep const indep1, xreg(2) coord(2) {p_end}
{phang2} matrix se46 = cov_nd {p_end}
{phang2} replace cutoff1 = 6 {p_end}
{phang2} replace cutoff2 = 4 {p_end}
{phang2} xgmlt C1 C2 cutoff1 cutoff2 newdep const indep1, xreg(2) coord(2) {p_end}
{phang2} matrix se64 = cov_nd {p_end}
{phang2} xgmlt C1 C2 cutoff1 cutoff2 newdep const indep1 instr2, xreg(3) coord(2) {p_end}

{title:Author}

{phang}
{cmd:Timbul Widodo} Olah Data Semarang.{break}
Homepage: {browse "http://www.youtube.com/@amalsedekah":http://www.youtube.com/@amalsedekah}. {break}
{p_end}