#!/usr/bin/env python3

"""
A simple wrapper of igv.js for visulization of BAM, BigWig and etc.
"""

import sys
import argparse
from os.path import basename
from os.path import isfile

def get_opt():
	group = argparse.ArgumentParser()

	group.add_argument("-r", "--ref", help="reference fasta file", required=True)
	group.add_argument("-m", "--bam", help="Input mapping file, in1.bam in2.bam ...", required=False, nargs = "*")
	group.add_argument("-w", "--bigwig", help="Input bigwig file, in1.bw in2.bw ...", required=False, nargs = "*")
	group.add_argument("-b", "--bed", help="bed annotation", required=False)
	group.add_argument("-g", "--gtf", help="gtf annotation", required=False, nargs = "*")
	
	return group.parse_args()

def build_bam_tracks(bams):
    tracks = []
    for bam in bams:
        if not isfile( bam + ".bai"):
            print("{}.bai is not existed".format(bam), file=sys.stderr)
            sys.exit()

        bam_id = basename(bam).replace(".bam","")
        track = """
        {{
            name: "{bam_id}",
            type: "alignment",
            format: "bam",
            url: "{bam}",
            indexURL: "{bam}.bai"
        }},""".format(bam = bam, bam_id = bam_id)
        tracks.append(track)

    tracks = "".join(tracks)[:-1]
    return tracks

def build_bw_tracks(bigwigs):
    tracks = []
    for bigwig in bigwigs:
        bw_id = basename(bigwig).replace(".bw","").replace(".bigwig","").replace(".bigWig","").replace("BigWig","")
        track = """
        {{
            name: "{bw_id}",
            format: "bigwig",
            url: "{bigwig}"
        }},""".format(bw_id = bw_id, bigwig = bigwig)
        tracks.append(track)

    tracks = "".join(tracks)[:-1]
    return tracks

def build_gtf_tracks(gtfs):
    tracks =  [] 
    for gtf in gtfs:
        track = """
	{{
	    type: "annotation",
            format: "gtf",
            sourceType: "file",
            url: "{gtf}",
            visibilityWindow: 500000,
            displayMode: "COLLAPSED",
            autoHeight: true
	}},""".format(gtf = gtf)
        tracks.append(track)
    tracks = "".join(tracks)[:-1]	
    return tracks

def build_bed_tracks(bed):
    track = """
	{{
	    type: "annotation",
            format: "bed",
            sourceType: "file",
            url: "{bed}",
            order: Number.MAX_VALUE,
            visibilityWindow: 300000000,
            displayMode: "EXPANDED"
	}}
	""".format(bed = bed)
    return track


def build_ref_track(fasta):
    if not isfile( fasta + ".fai"):
        print("{}.fai is not existed".format(fasta), file=sys.stderr)
        sys.exit()

    genome = basename(fasta)
    genome_id = genome.replace(".fasta","").replace(".fas","").replace(".fa","")

    bam_track ="""genome: "{genome}",
            reference: {{
                id: "{genome_id}",
                fastaURL: "{fasta}",
                indexURL: "{fasta}.fai"
            }}""".format(genome = genome, genome_id = genome_id, fasta = fasta)
    
    return bam_track


def make_html(genome, tracks):
    html = """<!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>igv.js</title>
   
        <script src="https://igv.org/web/release/2.6.6/dist/igv.min.js"></script>
    </head>

    <body>
        <div id="igv-div"></div>
    </body>
    <script>
        var igvDiv = document.getElementById("igv-div");
        var options = {{
            showNavigation: true,
            showRuler: true,
            {genome},
            tracks: [
                {tracks}
            ]
        }};
    
        igv.createBrowser(igvDiv, options)
            .then(function (browser) {{
                console.log("Created IGV browser");
            }})
    </script>
    </html>""".format(genome = genome, tracks = tracks)
    return html

def write_html(html):
    with open("index.html", "w") as file:
        file.writelines(html)

def igv_web(fasta, bams, bws, bed, gtfs):

    tracks = ""
    genome_track = build_ref_track(fasta)

    if bams is not None:
        bam_track = build_bam_tracks(bams)
        tracks += bam_track
    if bws is not None:
        bw_tracks = build_bw_tracks(bws)
        tracks = tracks + ",\n" + bw_tracks
    if gtfs is not None:
        gtf_track = build_gtf_tracks(gtfs)
        tracks = tracks + ",\n" + gtf_track
    if bed is not None:
        bed_track = build_bed_tracks(bed)
        tracks = tracks + ",\n" + bed_track

    html = make_html( genome_track, tracks)
    write_html(html)

if __name__ == "__main__":
     opts = get_opt()
     bams = opts.bam
     fasta = opts.ref
     bed = opts.bed
     gtfs = opts.gtf
     bws = opts.bigwig
     igv_web(fasta, bams, bws,  bed, gtfs)

