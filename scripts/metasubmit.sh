#!/bin/bash
#PBS -q hotel
#PBS -N metadensity
#PBS -l nodes=1:ppn=2
#PBS -l walltime=4:50:00
#PBS -o /home/hsher/cluster_msg/meta.out
#PBS -e /home/hsher/cluster_msg/meta.err
#PBS -t 0-218

source ~/miniconda3/etc/profile.d/conda.sh
conda activate metadensity

uids=(203 204 205 206 209 211 215 216 218 220 222 223 224 226 227 228 230 236 237 240 241 242 243 244 245 247 249 256 258 260 262 267 272 275 278 279 281 282 283 285 289 291 292 297 298 301 302 311 312 315 316 321 325 326 331 338 339 340 341 342 344 345 350 352 353 354 358 367 368 376 383 384 387 388 393 405 406 414 416 425 437 439 440 441 444 445 447 450 452 460 461 464 465 466 470 475 477 478 480 481 483 484 491 492 493 494 495 497 501 503 506 507 508 509 511 514 516 517 522 530 537 539 540 543 544 545 546 548 550 551 552 553 556 558 560 565 566 570 571 572 575 577 582 584 586 589 592 595 596 602 603 610 614 617 621 624 625 628 629 631 641 642 643 645 646 649 650 651 652 654 655 663 668 670 676 677 678 680 681 682 683 684 685 686 689 690 692 693 694 695 696 697 699 700 703 708 712 713 716 719 722reseq2 723reseq1 726 729 730 740 743 744 745 746 750 752 753 755 762 235x4000 284x4000fix 374x4000fix 632x)
uid=${uids[$PBS_ARRAYID]}

#python /home/hsher/projects/Metadensity/scripts/run_metadensity.py $uid ~/seqdata/metadensity

#python /home/hsher/projects/Metadensity/scripts/run_pos_enrich_encode3.py $uid 

#python /home/hsher/projects/Metadensity/scripts/run_kmer_from_read.py $uid 

#python /home/hsher/projects/Metadensity/scripts/run_shape_from_read.py $uid 

#python ~/Metadensity/scripts/rg4_enrichment.py $uid

python ~/Metadensity/scripts/run_metadensity.py $uid ~/densities/


