set terminal postscript eps color enhanced dashed size 16cm,12cm
set output "per_patient_plot.eps"
set multiplot layout 2,2 columnsfirst downwards
set ylabel "Relative Risk"
set xlabel "Patient"
#set xrange [0:15]
#set yrange [0:6]
set xtics ("1" 1, "2" 2, "3" 3, "4" 4, "5" 5, "6" 6, "7" 7, "8" 8, "9" 9, "10" 10, "Mean" 11.5, "R_{rbe}" 13, "R_n" 14)

set title "VMAT/C-ion bladder"
plot "-" using 1:2 with line lt 3 lw 1 lc rgb "#000000" notitle, \
     "-" using 1:2:3:4 with errorbars lt 1 lw 2 pt 13 ps 1.8 lc rgb "#1080FF" notitle, \
     "-" using 1:2:3:4 with errorbars lt 1 lw 2 pt 5 ps 1.4 lc rgb "#CE0000" notitle, \
     "-" using 1:2:3:4 with errorbars lt 1 lw 2 pt 6 ps 1.4 lc rgb "#5F5F5F" notitle
        0       1
        15      1
EOF
        1       1.171645778     0.961413018     1.394799109
        2       0.604246864     0.312390561     1.248998076
        3       0.988786061     0.724358007     1.418965072
        4       1.624472965     1.438678692     1.826082772
        5       0.714958052     0.415706964     1.254621617
        6       1.198187562     0.946193671     1.668772585
        7       0.333572754     0.101592780     1.021343314
        8       1.063671254     0.746502158     1.677268075
        9       4.707986164     3.293797052     5.655876101
        10      0.711300908     0.466789412     1.144087787

EOF
        11.5    1.3118837       0.6505619       2.1848833
EOF
        13      1.19138122      0.742120694     1.60496705
        14      0.818930211     0.434227531     1.2823225
EOF

set title "VMAT/C-ion rectum"
plot "-" using 1:2 with line lt 3 lw 1 lc rgb "#000000" notitle, \
     "-" using 1:2:3:4 with errorbars lt 1 lw 2 pt 13 ps 1.8 lc rgb "#1080FF" notitle, \
     "-" using 1:2:3:4 with errorbars lt 1 lw 2 pt 5 ps 1.4 lc rgb "#CE0000" notitle, \
     "-" using 1:2:3:4 with errorbars lt 1 lw 2 pt 6 ps 1.4 lc rgb "#5F5F5F" notitle
        0       1
        15      1
EOF
        1       0.731211612     0.704569515     0.842373933
        2       0.755719433     0.568080745     1.118121185
        3       0.748155656     0.592323138     1.028420393
        4       0.341391034     0.281488683     0.493427448
        5       0.752358919     0.627301135     0.859249171
        6       0.475063482     0.410892189     0.649106292
        7       0.136445525     0.068345431     0.282798611
        8       0.593259324     0.407947289     0.975271154
        9       0.700134932     0.661364996     0.734153826
        10      0.584453383     0.479882228     0.754444698
EOF
        11.5    0.5817372       0.4129403       0.7991492
EOF
        13      0.691522274     0.456031126     1.00349199
        14      0.433118198     0.299480126     0.572087099
EOF

set title "VMAT/IMPT bladder"
plot "-" using 1:2 with line lt 3 lw 1 lc rgb "#000000" notitle, \
     "-" using 1:2:3:4 with errorbars lt 1 lw 2 pt 13 ps 1.8 lc rgb "#1080FF" notitle, \
     "-" using 1:2:3:4 with errorbars lt 1 lw 2 pt 5 ps 1.4 lc rgb "#CE0000" notitle, \
     "-" using 1:2:3:4 with errorbars lt 1 lw 2 pt 6 ps 1.4 lc rgb "#5F5F5F" notitle
        0       1
        15      1
EOF
        1       2.183086905     1.812047142     2.377832054
        2       0.974569024     0.518517680     1.820262433
        3       1.972346324     1.498271617     2.416930492
        4       3.001395504     2.649911528     3.331453128
        5       1.172478488     0.685117361     1.853291487
        6       2.817014904     2.256998818     3.331149334
        7       0.492619711     0.155769388     1.385824767
        8       2.309263426     1.720686034     2.982674746
        9       1.222529947     0.864555791     1.423381843
        10      1.047480551     0.693279665     1.505378340
EOF
        11.5    1.7191894       1.0638891       2.3668484
EOF
        13      1.70656419      1.0524295       2.35775285
        14      1.44391866      0.762823452     2.22061246
EOF

set title "VMAT/IMPT rectum"
plot "-" using 1:2 with line lt 3 lw 1 lc rgb "#000000" notitle, \
     "-" using 1:2:3:4 with errorbars lt 1 lw 2 pt 13 ps 1.8 lc rgb "#1080FF" notitle, \
     "-" using 1:2:3:4 with errorbars lt 1 lw 2 pt 5 ps 1.4 lc rgb "#CE0000" notitle, \
     "-" using 1:2:3:4 with errorbars lt 1 lw 2 pt 6 ps 1.4 lc rgb "#5F5F5F" notitle
        0       1
        15      1
EOF
        1       1.602046599     1.496764648     1.733400791
        2       1.667402242     1.364943901     2.020864632
        3       1.565453852     1.321213085     1.804709332
        4       0.575080612     0.492440283     0.743089466
        5       1.307116713     1.142305598     1.377834874
        6       0.865952338     0.782636229     1.027425695
        7       0.182971973     0.094254371     0.359495367
        8       1.282696740     0.953735923     1.782552448
        9       0.953971999     0.883230427     0.992639503
        10      0.974469762     0.834651001     1.116190719
EOF
        11.5    1.0981823       0.7773345       1.4253127
EOF
        13      1.0745199       0.757330632     1.40408845
        14      1.0305244       0.676997974     1.40034721
EOF
