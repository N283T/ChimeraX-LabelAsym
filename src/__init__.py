"""ChimeraX-LabelAsym: expose mmCIF label_asym_id for residue selection.

After this bundle initializes, each residue of any structure opened from
an mmCIF file (or mmCIF-derived metadata) gets a ``label_asym_id`` string
attribute. Users can then select via::

    select ::label_asym_id=A
    color ::label_asym_id=B red
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


bundle_api = _LabelAsymAPI()
