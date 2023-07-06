#!/bin/bash

rclone copy --progress --transfers 12 --checkers 12 --copy-links ./results/ gemini32:/mnt/gemini3-4/PROJECTS/new_RNA_isoform_expression_across_tissues
