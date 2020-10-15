#!/bin/csh
#PBS -q hotel
#PBS -N metadensity
#PBS -l nodes=1:ppn=6
#PBS -l walltime=4:50:00
#PBS -o /home/hsher/cluster_msg/meta.out
#PBS -e /home/hsher/cluster_msg/meta.err
#PBS -t 0-102

set uids=(203 204 205 206 209 211 215 216 218 222 223 227 228 230 272 281 282 283 291 292 297 298 301 302 311 321 344 345 383 384 387 388 393 437 477 478 484 492 493 494 495 497 501 503 539 540 543 544 545 546 552 553 556 558 571 572 575 589 592 593 595 596 602 603 617 621 625 628 629 631 641 642 643 649 650 651 652 654 655 678 692 694 695 696 703 708 712 713 722reseq2 723reseq1 726 729 730 743 744 745 750 752 753 755 762 235x4000 632x)
set uid=$uids[$PBS_ARRAYID]

python /home/hsher/projects/Metadensity/scripts/run_metadensity.py $uid ~/seqdata/metadensity


