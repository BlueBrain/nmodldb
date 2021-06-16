class MechInfo:
    def __init__(self, name):
        self.name = name
        self.nlines = 0
        self.nverbatim = 0


class ChannelInfo(MechInfo):
    def __init__(self, name):
        MechInfo.__init__(self, name)
        self.is_channel = True


class SynapseInfo(MechInfo):
    def __init__(self, name):
        MechInfo.__init__(self, name)
        self.is_synapse = True


class ModelInfo:
    def __init__(self, url):
        self.name = None
        self.description = ""
        self.citation = ""
        self.neuron_model = False
        self.simulator = None
        self.url = url

        self.is_zip = False
        self.is_git = False
        self.git_branch = None

        self.output_dir = None
        self.cache_dir = None
        self.output_filename = None

        self.mod_dir = "."
        self.mod_filter = None

        self.nmod = 0
        self.nchannels = 0
        self.nsynapses = 0
        self.nverbatim = 0

        self.channels = []
        self.synapses = []

        self.group = None
