import click
from flask_admin.contrib.sqla import ModelView

from .tomo_table_def import TomoConfig


### COMMAND LINE INTERFACE PLUGIN DEFINITION

@click.group()
def tomo():
    """Package to compute dispersion curves using automated FTAN and
    invert them to obtain a velocity model at different periods."""
    pass

@click.command()
def info():
    from msnoise.api import connect, get_config
    from .default import default
    db = connect()
    click.echo('')
    click.echo('Raw config bits: "D"efault or "M"odified (green)')
    for key in default.keys():
        tmp = get_config(db, key, plugin='Tomo')
        if tmp == default[key][1]:
            click.secho(" D %s: %s" %(key, tmp ))
        else:
            click.secho(" M %s: %s" %(key, tmp ), fg='green')

@click.command()
def install():
    from .install import main
    main()


@click.command()
def prepare_ccf():
    from .export_single_sided import main
    main()


@click.command()
def ftan_example():
    from .examplepickdispcurve import main
    main()

@click.option('-p', '--pair', default=None,  help='FTAN a specific pair',
              multiple=True)
@click.option('-b', '--bmin', default=None,  help='force bmin',)
@click.option('-B', '--bmax', default=None,  help='force bmax',)
@click.option('-s', '--show', default=1,  help='show plot',)
@click.command()
def ftan(pair, bmin, bmax, show):
    from .ftan import main
    main(pair, bmin, bmax, show)


@click.command()
def prepare_tomo():
    from .prepare_tomo import main
    main()

@click.option('-p', '--per',type=float, default=None,  help='force per',)
@click.option('--a1', type=float, default=None, help='force bmin',)
@click.option('--b1', type=float, default=None, help='force bmin',)
@click.option('--l1', type=float, default=None, help='force bmin',)
@click.option('--s1', type=float, default=None, help='force bmin',)
@click.option('--a2', type=float, default=None, help='force bmin',)
@click.option('--b2', type=float, default=None, help='force bmin',)
@click.option('--l2', type=float, default=None, help='force bmin',)
@click.option('--s2', type=float, default=None, help='force bmin',)
@click.option('-s', '--show', help='Show interactively?',
              default=True, type=bool)
@click.command()
def answt(per, a1, b1, l1, s1, a2, b2, l2, s2, show):
    from .ANSWT import main
    main(per, a1, b1, l1, s1, a2, b2, l2, s2, show)


@click.command()
def prepare_1d():
    from .prepare_1d import main
    main()


@click.command()
def plot():
    from .plotdisp import main
    main()

@click.command()
def plot3d():
    from .plot3d import main
    main()


tomo.add_command(info)
tomo.add_command(ftan_example)
tomo.add_command(prepare_ccf)
tomo.add_command(prepare_tomo)
tomo.add_command(ftan)
tomo.add_command(install)
tomo.add_command(answt)
tomo.add_command(plot)
tomo.add_command(prepare_1d)
tomo.add_command(plot3d)


### WEB INTERFACE PLUGIN DEFINITION
class TomoConfigView(ModelView):
    # Disable model creation
    view_title = "MSNoise TOMO Configuration"
    name = "Configuration"

    #inline_models = (SaraConfig,)
    can_create = False
    can_delete = False
    page_size = 50
    # Override displayed fields
    column_list = ('name', 'value')

    def __init__(self, session, **kwargs):
        # You can pass name and other parameters if you want to
        super(TomoConfigView, self).__init__(TomoConfig, session,
                                             endpoint="tomoconfig",
                                             name="Config",
                                             category="Tomo", **kwargs)

def getitem(obj, item, default):
    if item not in obj:
        return default
    else:
        return obj[item]


# Job definitions

def register_job_types():
    jobtypes = []
    jobtypes.append( {"name":"TOMO_SAC", "after":"refstack"} )
    return jobtypes
