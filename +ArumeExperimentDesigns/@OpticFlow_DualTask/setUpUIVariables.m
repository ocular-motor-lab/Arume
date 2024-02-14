function this = setUpUIVariables(this)

    % fixation rect stuff
    hfixwidth = 3;
    hfixlen = 40;
    fix1 = [this.Graph.pxWidth/2-hfixwidth;
        this.Graph.pxHeight/2-hfixlen;
        this.Graph.pxWidth/2+hfixwidth;
        this.Graph.pxHeight/2+hfixlen];
    fix2 =  [this.Graph.pxWidth/2-hfixlen;
        this.Graph.pxHeight/2-hfixwidth;
        this.Graph.pxWidth/2+hfixlen;
        this.Graph.pxHeight/2+hfixwidth];
    this.uicomponents.fixrects = [fix1,fix2];
    
    % initialize the anchor point for the response line
    this.uicomponents.resplinex1 = this.Graph.pxWidth/2;
    this.uicomponents.respliney1 =this.Graph.pxHeight/2;
    this.uicomponents.linelen = this.Graph.pxHeight/8;
    this.uicomponents.anginc = .01;
    
    % angle ticks
    rulerlineangs = linspace(0,pi,5);
    rulerlinex1 = this.uicomponents.resplinex1+cos(rulerlineangs)*this.uicomponents.linelen*.9;
    rulerlinex2 = this.uicomponents.resplinex1+cos(rulerlineangs)*this.uicomponents.linelen*1;
    rulerliney1 = this.uicomponents.respliney1-sin(rulerlineangs)*this.uicomponents.linelen*.9;
    rulerliney2 = this.uicomponents.respliney1-sin(rulerlineangs)*this.uicomponents.linelen*1;
    this.uicomponents.rulerlinexy = [reshape([rulerlinex1;rulerlinex2],1,[]);...
        reshape([rulerliney1;rulerliney2],1,[])];
    
    % outside semi-circle
    outsidelineangs = linspace(0,pi,100);
    outsidelinex = this.uicomponents.resplinex1+cos(outsidelineangs)*this.uicomponents.linelen;
    outsidelinex = [outsidelinex(1),repelem(outsidelinex(2:end),2),outsidelinex(1)];
    
    outsideliney = this.uicomponents.respliney1-sin(outsidelineangs)*this.uicomponents.linelen;
    outsideliney = [outsideliney(1),repelem(outsideliney(2:end),2),outsideliney(1)];
    
    this.uicomponents.outsidelinexy = [outsidelinex;outsideliney];
    
    % location for degree text box response
    this.uicomponents.texty = this.Graph.pxHeight/2+this.Graph.pxHeight/20;

    % text box response for visual search bounting task
    % this.uicomponents.visualsearchtexty1 = this.Graph.pxHeight/2-this.Graph.pxHeight/20;
    % this.uicomponents.visualsearchtexty2 = this.Graph.pxHeight/2-this.Graph.pxHeight/20;

end