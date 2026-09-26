#!/usr/bin/env pvpython
"""
Animate the tracer particles of a run around its (static) obstacle in ParaView.

The field file the solver writes at the end of a run carries `fluidFraction`,
0 in solid and 1 in open fluid; its 0.5 contour is the body. The particle files
in vis_files/ are indexed by vis_files/particles.vtk.series, which gives every
file its simulation time. Both are loaded here and set up as one scene.

Run from the directory the solver ran in:

    paraview --script=tools/paraview/animate-particles.py
    pvbatch  tools/paraview/animate-particles.py karman.vtk --save karman.avi
    pvbatch  tools/paraview/animate-particles.py karman.vtk --save frame.png

The first opens the GUI with the scene ready to press Play; so does View >
Python Shell > Run Script in a GUI already open. Neither passes arguments, so
the field file then defaults to the only *.vtk in the working directory.

--save with a movie extension (.avi, .ogv) or an image extension (.png, .jpg)
writes the whole animation; an image extension writes one numbered frame per
time step. --frame picks a single time instead.
"""
import argparse
import glob
import os
import sys

from paraview.simple import (
    ColorBy,
    Contour,
    GetActiveViewOrCreate,
    GetAnimationScene,
    GetTimeKeeper,
    LegacyVTKReader,
    OpenDataFile,
    Outline,
    Render,
    ResetCamera,
    SaveAnimation,
    SaveScreenshot,
    Show,
)


def parse(argv):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("field", nargs="?", help="field output of the run, e.g. karman.vtk")
    ap.add_argument("--particles", default="vis_files",
                    help="directory holding the particle files (default vis_files)")
    ap.add_argument("--save", help="write the animation (.avi, .ogv, .png, ...)")
    ap.add_argument("--frame", type=float,
                    help="with --save, write only the frame at this time")
    ap.add_argument("--size", type=int, nargs=2, default=[1600, 800],
                    metavar=("W", "H"), help="image size (default 1600 800)")
    ap.add_argument("--point-size", type=float, default=3.0,
                    help="particle point size in pixels (default 3)")
    # Inside the GUI the arguments are ParaView's own, so ignore what is not ours.
    args, _ = ap.parse_known_args(argv)
    if args.field and not (args.field.endswith(".vtk") and os.path.exists(args.field)):
        args.field = None
    return args


def field_file(given):
    if given:
        return given
    found = sorted(glob.glob("*.vtk"))
    if len(found) != 1:
        sys.exit("give the field file: found %d *.vtk here" % len(found))
    return found[0]


def particle_source(directory):
    series = os.path.join(directory, "particles.vtk.series")
    if os.path.exists(series):
        return OpenDataFile(series)

    # Written by an older solver: the files still animate, by file index.
    files = sorted(glob.glob(os.path.join(directory, "particles_*.vtk")))
    if not files:
        sys.exit("no particle files in %s" % directory)
    print("no %s, animating by file index" % series)
    return LegacyVTKReader(FileNames=files)


def main(argv):
    args = parse(argv)
    view = GetActiveViewOrCreate("RenderView")
    view.ViewSize = args.size
    view.OrientationAxesVisibility = 1
    view.Background = [1.0, 1.0, 1.0]
    try:
        # ParaView 5.10 and later pick the background from a palette unless told.
        view.UseColorPaletteForBackground = 0
    except AttributeError:
        pass

    field = LegacyVTKReader(FileNames=[field_file(args.field)])
    field.UpdatePipeline()

    outline = Show(Outline(Input=field), view)
    outline.DiffuseColor = [0.2, 0.2, 0.2]
    outline.AmbientColor = [0.2, 0.2, 0.2]

    names = field.PointData.keys()
    if "fluidFraction" in names:
        body = Contour(Input=field, ContourBy=["POINTS", "fluidFraction"],
                       Isosurfaces=[0.5], ComputeScalars=0)
        body.UpdatePipeline()
        if body.GetDataInformation().GetNumberOfPoints() > 0:
            shown = Show(body, view)
            ColorBy(shown, None)
            shown.DiffuseColor = [0.6, 0.6, 0.65]
            shown.AmbientColor = [0.6, 0.6, 0.65]
            shown.Specular = 0.3
        else:
            print("fluidFraction has no 0.5 contour: the domain has no obstacle")
    else:
        print("%s has no fluidFraction, so no obstacle is drawn" % args.field)

    particles = particle_source(args.particles)
    dots = Show(particles, view)
    dots.SetRepresentationType("Points")
    dots.PointSize = args.point_size
    dots.RenderPointsAsSpheres = 1
    dots.DiffuseColor = [0.1, 0.35, 0.8]
    dots.AmbientColor = [0.1, 0.35, 0.8]

    scene = GetAnimationScene()
    scene.UpdateAnimationUsingDataTimeSteps()
    times = GetTimeKeeper().TimestepValues

    # Looking obliquely down on the flow from upstream, the usual view of a wake.
    ResetCamera(view)
    camera = view.GetActiveCamera()
    camera.Elevation(25)
    camera.Azimuth(-30)
    view.ResetCamera()
    camera.Dolly(1.3)

    if times:
        scene.AnimationTime = times[0]
    Render(view)

    if args.save:
        if args.frame is not None:
            scene.AnimationTime = args.frame
            Render(view)
            SaveScreenshot(args.save, view, ImageResolution=args.size)
        else:
            SaveAnimation(args.save, view, ImageResolution=args.size)
        print("wrote %s" % args.save)


# Not guarded by __name__: the GUI runs a script under a name of its own.
main(sys.argv[1:])
