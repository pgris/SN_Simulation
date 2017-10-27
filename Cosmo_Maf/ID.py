class Id:
    def __init__(self,thedir,thedir_obs,fieldname,fieldid,X1,Color,season,colorfig):
        self.thedir=thedir
        self.thedir_obs=thedir_obs
        self.fieldname=fieldname
        self.fieldid=fieldid
        self.X1=X1
        self.Color=Color
        self.season=season
        self.colorfig=colorfig

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
    def colorfig(self):
        return self.colorfig
