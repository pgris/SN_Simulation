class Id:
    def __init__(self,thedir,thedir_obs,fieldname,fieldid,X1,Color,season,z,colorfig,lstyle,T0step,simul_name):
        self.thedir=thedir
        self.thedir_obs=thedir_obs
        self.fieldname=fieldname
        self.fieldid=fieldid
        self.X1=X1
        self.Color=Color
        self.season=season
        self.colorfig=colorfig
        self.z=z
        self.linestyle=lstyle
        self.T0step=T0step
        self.simul_name=simul_name

    @property
    def thedir(self):
        return self.thedir
    @property
    def thedir_obs(self):
        return self.thedir_obs
    @property
    def fieldname(self):
        return self.fieldname
    @property
    def fieldid(self):
        return self.fieldid
    @property
    def X1(self):
        return self.X1
    @property
    def Color(self):
        return self.Color
    @property
    def season(self):
        return self.season
    @property
    def z(self):
        return self.z
    @property
    def colorfig(self):
        return self.colorfig
    @property
    def linestyle(self):
        return self.linestyle
    @property
    def T0step(self):
        return self.T0step
    @property
    def simul_name(self):
        return self.simul_name
