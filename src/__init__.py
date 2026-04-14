"""ChimeraX-LabelAsym: expose mmCIF label_asym_id for residue selection.

After this bundle initializes, each residue of any structure opened from
an mmCIF file (or mmCIF-derived metadata) gets a ``label_asym_id`` string
attribute, plus dynamic selectors ``la_<label>``. Users can select via::

    select ::label_asym_id="A"   # attribute form (any label)
    select la_E                   # short form (registered on open)
    labelcolor                    # color atoms by label_asym_id
"""

from chimerax.core.toolshed import BundleAPI


class _LabelAsymAPI(BundleAPI):
    api_version = 1

    @staticmethod
    def initialize(session, bundle_info):
        from .hook import install

        install(session)

    @staticmethod
    def finish(session, bundle_info):
        from .hook import uninstall

        uninstall(session)

    @staticmethod
    def register_command(bundle_info, command_info, logger):
        from .commands import register_commands

        register_commands(logger)


bundle_api = _LabelAsymAPI()
