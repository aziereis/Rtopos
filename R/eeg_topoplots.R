#' A topoplot Function
#'
#' This function allows you to plot ERP topoplots for different conditions, specify the printed limits (microvolt scale) and customize the appearance (colours, fonts, etc.).
#' @param dataset the data set to plot
#' @param condition name of the variable that specifies the condition to plot
#' @param conditionlevels if not specified, all levels available from the unique values of the condition will be used. (currently max 4)
#' @param differenceplots calculate differences between condition levels and plot them. Defaults to TRUE.
#' @param differencelevel specify a specific condition as the reference level. Defaults to "all" for all pairwise comparisons.
#' @param extralegend plot legend extra (vs. under the label) to the figure. Defaults to TRUE.
#' @param addlabels include the labels of the electrode channels. Defaults to FALSE.
#' @param method one of "akima" or "biharmonic". Akima performs faster and includes a spline method but requires a specific version of the akima package (not newer than 0.6-2.1). Defaults to "akima".
#' @param simsignals Defaults to TRUE.
#' @param quick faster plot with lower resolution. Defaults to TRUE.
#' @param palette colour palette of the topoplot (e.g., BrBG, PiYG, PRGn, PuOr, RdBu, RdGy, RdYlBu, RdYlGn, Spectral). Defaults to "RdBu".
#' @param signallimits lmit the colour palette signals for all plots to make them comparable. Defaults to NULL.
#' @param signalcomplimits same as signallimits but for the difference plots. Defaults to NULL.
#' @param markchans a character vector specifying the names of the to be highlighted ROI channels. Defaults to NULL.
#' @param nrchans choose 64 or 128. 64 includes the extended 10-20 system layour, 128 the ABCD layout. Defaults to 64.
#' @param chanmarkcol colour of the ROIchannels if to highlight. Defaults to "black".
#' @param plottitle title of the plot, Defaults to ''
#' @param marksize size of the electrode markers. Defaults to 1.1
#' @param textsize size of the text capture
#' @param textfamily if specified, the package "extrafont" needs to be installed and loaded.
#' @keywords topo
#' @export
#' @examples
#'
#'
#'
#'
#' data(xdat) #load example data frame
#' maketopoplot(dataset = xdat, condition = "difficulty", nrchans = 128)
maketopoplot<-function(
    dataset=dat,
    condition=NULL,
    conditionlevels=levels(as.factor(dataset[,condition])),
    differenceplots=T,
    differencelevel = "all",
    extralegend = T,
    addlabels=F,
    method = "akima",
    simsignals=T,
    quick=T,
    palette = "RdBu",
    signallimits=NULL,
    signalcomplimits =NULL,
    markchans = NULL,
    nrchans = 64,
    chanmarkcol = "black",
    plottitle='',
    marksize = 1.1,
    textsize = 10,
    textfamily = NULL){
  #Diverging
  #BrBG, PiYG, PRGn, PuOr, RdBu, RdGy, RdYlBu, RdYlGn, Spectral

  ##todebug
  # dataset=as.data.frame(xdat)
  # condition="difficulty"
  # conditionlevels=levels(as.factor(dataset[,condition]))
  # differenceplots = T
  # differencelevel = "all"
  # extralegend = T
  # palette="custom"
  # addlabels = F
  # quick=T
  # signallimits = NULL
  # simsignals=T
  # plottitle=''
  # nrchans = 128
  # palette = "RdBu"
  # signalcomplimits =NULL
  # markchans = NULL
  # chanmarkcol = "black"
  # textsize = 10
  # method = "akima"
  # textfamily ="Arial"
  #load packages that are needed

  akimaversion = as.character(packageVersion("akima"))
  if(akimaversion !="0.6.2.1"){
    method = "biharmonic"
  }

  limsignals = T
  dataset<-as.data.frame(dataset)
  #check conditions to plot
  if(length(condition)==0){
    condition = "nocondition"
    dataset$nocondition<-factor(1)
  }

  if (differencelevel != "all"){
    if(differencelevel%in%conditionlevels){
      dataset[,c(condition)]<-relevel(as.factor(dataset[,c(condition)]), ref = paste0(differencelevel))
      levelidx<-match(differencelevel, conditionlevels)
      conditionlevels = c(conditionlevels[levelidx], conditionlevels[-c(levelidx)])
    }
  }
  dataset[c(condition)]<-lapply(dataset[c(condition)], factor)

  if(class(conditionlevels)=="NULL"){message("Warning. something went wrong with conditions. check whether the dataset has the correct format. conditions as factors")}


  originalconditionlevels<-conditionlevels #store the original conditions (later there will be new ones added (e.g. comparisonplots))


  abbrvlabel<-function(x){
    x = as.character(x)
    if (nchar(x) <= 7){
      return(x)
    }
    else{
      return(abbreviate(x,7, method ="both.sides"))
    }
  }

  #circle fun to make the head and to exclude values outside the head
  circleFun <- function(center = c(0,0),diameter = 1, npoints = 100){
    r = diameter / 2
    tt <- seq(0,2*pi,length.out = npoints)
    xx <- center[1] + r * cos(tt)
    yy <- center[2] + r * sin(tt)
    return(data.frame(x = xx, y = yy))
  }
  circledat <- circleFun(c(0.5, 0.5), 1, npoints = 10000) # center on [.5, .5] #black circle

  if(quick){
    circledat2 <- circleFun(c(0.5, 0.5), 1.02, npoints = 10000) # center on [.5, .5] #white circle for fast plots (cut points away)
    circledat3 <- circleFun(c(0.5, 0.5), 1.03, npoints = 10000) # center on [.5, .5] #white circle for fast plots (cut points away)
  }
  if(!quick){
    circledat2 <- circleFun(c(0.5, 0.5), 1.03, npoints = 10000) # center on [.5, .5] #white circle for fast plots (cut points away)
  }



  #Load electrode locations
  #Load electrode locations
  if (nrchans == 64){
    #electrodelocs<-read.delim("ElectrodeLocation64_fieldtrip.txt", header=F, sep=" ")
    electrodelocs<-Rtopos::electrodelocs64
    electrodelocs[,c(1,5)]<-NULL
    names(electrodelocs)<-c("x", "y", "z","channel")
    electrodelocs$channel<-tolower(electrodelocs$channel)
    electrodelocs[,1:3]<-electrodelocs[,1:3]+0.5
    # if I should only plot one condition and one level of that condition (e.g. emotion - happy)
    # then it doesn't make sense to do any difference plots
    if(length(conditionlevels)==1){
      differenceplots=F}
    #now load and aggregate data
    electrodechanlist<-electrodelocs$channel
    dat<-tidyr::pivot_longer(dataset, cols=matches(electrodechanlist), names_to="channel", values_to="signal") #convert to long format
  }
  if (nrchans == 128){

    electrodelocs<-Rtopos::electrodelocs128#read.delim("ElectrodeLocation128_fieldtrip.txt", header=F, sep=" ")
    electrodelocs[,c(1,5)]<-NULL
    names(electrodelocs)<-c("x", "y", "z","channel")
    electrodelocs$channel<-tolower(electrodelocs$channel)
    electrodelocs[,1:3]<-electrodelocs[,1:3]+0.5


    # if I should only plot one condition and one level of that condition (e.g. emotion - happy)
    # then it doesn't make sense to do any difference plots
    if(length(conditionlevels)==1){
      differenceplots=F}

    #now load and aggregate data
    electrodechanlist<-electrodelocs$channel
    dat<-tidyr::pivot_longer(dataset, cols=matches(electrodechanlist), names_to="channel", values_to="signal") #convert to long format

  }


  # if I should only plot one condition and one level of that condition (e.g. emotion - happy)
  # then it doesn't make sense to do any difference plots
  if(length(conditionlevels)==1){
    differenceplots=F}


  dat<-as.data.frame(dat) #make a dataframe out of it (there were some compatibility problems if it stays a tibble)

  #only select those conditions to analyse
  dat<-subset(dat, dat[,condition] %in% conditionlevels )

  dat2<-dat%>%dplyr::group_by(channel, !!rlang::sym(condition))%>%dplyr::summarize(signal=mean(signal, na.rm=T)) #calculate the mean signal for each condition and channel
  #length(conditionlevels) #check how many conditionlevels to plot

  #if there are 2 conditionlevels and I should do difference plots
  if(length(conditionlevels)==2 & differenceplots==T){
    dat22<-tidyr::spread(dat2, !!rlang::sym(condition), signal)
    dat22<-dat22%>%dplyr::mutate(difference = get(conditionlevels[2])-get(conditionlevels[1]))
  }

  #if there are 3 conditionlevels and I should do difference plots
  if(length(conditionlevels)==3 & differenceplots==T){
    dat22<-tidyr::spread(dat2, !!rlang::sym(condition), signal)
    if (differencelevel =="all"){
      dat22<-dat22%>%dplyr::mutate(comp1 = get(conditionlevels[2])-get(conditionlevels[1]),
                                   comp2 = get(conditionlevels[3])-get(conditionlevels[1]),
                                   comp3 = get(conditionlevels[2])-get(conditionlevels[3]))
    }else{
      dat22<-dat22%>%dplyr::mutate(comp1 = get(conditionlevels[2])-get(conditionlevels[1]),
                                   comp2 = get(conditionlevels[3])-get(conditionlevels[1]))
    }


  }
  if(length(conditionlevels)==4 & differenceplots==T){
    differencelevel=conditionlevels[1]
    message("not all pairwise differences will be plotted")
    dat22<-tidyr::spread(dat2, !!rlang::sym(condition), signal)
    dat22<-dat22%>%dplyr::mutate(comp1 = get(conditionlevels[2])-get(conditionlevels[1]),
                                 comp2 = get(conditionlevels[3])-get(conditionlevels[1]),
                                 comp3 = get(conditionlevels[4])-get(conditionlevels[1]))



  }
  if(length(conditionlevels)==5 & differenceplots==T){
    differencelevel=conditionlevels[1]
    message("not all pairwise differences will be plotted")
    dat22<-tidyr::spread(dat2, !!rlang::sym(condition), signal)
    dat22<-dat22%>%dplyr::mutate(comp1 = get(conditionlevels[2])-get(conditionlevels[1]),
                                 comp2 = get(conditionlevels[3])-get(conditionlevels[1]),
                                 comp3 = get(conditionlevels[4])-get(conditionlevels[1]),
                                 comp4 = get(conditionlevels[5])-get(conditionlevels[1]))



  }
  if(length(conditionlevels)==6 & differenceplots==T){
    message("not all pairwise differences will be plotted")
    dat22<-tidyr::spread(dat2, !!rlang::sym(condition), signal)
    dat22<-dat22%>%dplyr::mutate(comp1 = get(conditionlevels[2])-get(conditionlevels[1]),
                                 comp2 = get(conditionlevels[3])-get(conditionlevels[1]),
                                 comp3 = get(conditionlevels[4])-get(conditionlevels[1]),
                                 comp4 = get(conditionlevels[5])-get(conditionlevels[1]),
                                 comp5 = get(conditionlevels[6])-get(conditionlevels[1]))

  }
  if(length(conditionlevels)>1 & differenceplots==T){
    dat2<-tidyr::gather(dat22, condition, signal,  -channel, factor_key=TRUE)
    if(length(conditionlevels)>2 & differenceplots==T){
      if(differencelevel=="all"){
        comp1 = paste0(abbrvlabel(conditionlevels[2])," - ",abbrvlabel(conditionlevels[1]))
        comp2 = paste0(abbrvlabel(conditionlevels[3])," - ",abbrvlabel(conditionlevels[1]))
        comp3 = paste0(abbrvlabel(conditionlevels[2])," - ",abbrvlabel(conditionlevels[3]))
      }else{
        for(idx in 2:length(originalconditionlevels)){
          assign(paste0("comp",idx-1), paste0(abbrvlabel(originalconditionlevels[idx])," - ",abbrvlabel(originalconditionlevels[1])))

        }
      }
    }
  }


  #prepare dataset
  names(dat2)<- c("channel","condition", "signal")
  dat2<-droplevels(dat2) # drop unused factor levels- important. otherwise there will be errors

  (conditionlevels<-unique(dat2$condition)) # to also include the difference plots)

  # now calculate the topoplots based on how many levels to compare and plot
  plotlist<-list()  #create empty list for plots
  legendlist<-list() #create empty list seperately for legends
  signalss<-list()  #create empty list for minimum and maximum signals (needed for limit signals)

  interpdataset<-list() #create empty list with the data that can later be put back into the loop
  levelz = conditionlevels[1]
  for(levelz in conditionlevels){

    message(levelz) #which level is to be calculated
    dat<-subset(dat2, condition==levelz, select=c(channel, signal))
    dat<-merge(dat, electrodelocs, by="channel")
    #class(dat) # just to check

    #length(dat[is.na(dat)])
    dat<-na.omit(dat)
    #length(dat[is.na(dat)])
    #should the plot be run more quick but the edges are not so nicely defined?
    if(quick){
      points2int<-100
    }else{
      points2int<-500
    }
    #here the interpolation and extrapolation happens:
    if(method =="akima"){
    datmat <- akima::interp(dat$x, dat$y, dat$signal,
                            xo = seq(0, 1, length = points2int),
                            yo = seq(0, 1, length = points2int),
                            extrap=TRUE, linear=FALSE, duplicate = "error")
    datmat2 <- data.frame(
      do.call(expand.grid, lapply(dim(datmat$z), seq_len)),
      value = as.vector(datmat$z)
    )
    }

    if(method =="biharmonic"){
      datmat2 <-Rtopos:::biharmonic(dat,
                           grid_res = points2int)
    }

    names(datmat2) <- c('x', 'y', 'value')
    datmat2[,1:2] <- datmat2[,1:2]/points2int
    # ignore anything outside the circle
    datmat2$incircle <- (datmat2$x - .5)^2 + (datmat2$y - .5)^2 < .5^2 # mark

    datmat2 <- datmat2[datmat2$incircle,]

    interpdataset[[levelz]]<-datmat2 #store data in the list
    signalss[[levelz]]<-c(min(datmat2$value),max(datmat2$value)) # retrieve minimum and maximum signals (inside the circle)
  }



  #signalrange<-unlist(signalss) #for all conditions
  comparisons<-base::setdiff(conditionlevels,originalconditionlevels) # which conditionlevels are differences?
  #if signals should be limited and there are difference plots
  if(limsignals==T & differenceplots==T & class(signallimits)=="NULL"){
    minsignalcomp<-min(unlist(signalss[names(signalss)%in%comparisons]), na.rm =T)#for the difference plots
    maxsignalcomp<-max(unlist(signalss[names(signalss)%in%comparisons]), na.rm =T)
  }

  # loop through conditions ----
  #get the signal limits (either calculate them based on data or take those that are provided by the user)
  for(levelz in conditionlevels){
    if(limsignals==T & levelz%in%originalconditionlevels & class(signallimits)=="NULL"){
      minsignal<-min(unlist(signalss[names(signalss)%in%originalconditionlevels]))
      maxsignal<-max(unlist(signalss[names(signalss)%in%originalconditionlevels]))
      if(simsignals ==T){
        maxabssignal<-max(abs(c(minsignal, maxsignal)))
        minsignal<--(maxabssignal)
        maxsignal<-maxabssignal
      }
    }
    if(limsignals==T & levelz%in%comparisons & class(signallimits)=="NULL"){
      minsignal<-minsignalcomp
      maxsignal<-maxsignalcomp
      if(simsignals ==T){
        maxabssignal<-max(abs(c(minsignal, maxsignal)))
        minsignal<--(maxabssignal)
        maxsignal<-maxabssignal
      }
      if(class(signalcomplimits)!="NULL"){
        minsignal<-min(signalcomplimits)
        maxsignal<-max(signalcomplimits)
      }
    }

    #user provided limits:
    if(class(signallimits)!="NULL"){
      if(limsignals==T & levelz%in%originalconditionlevels){
        minsignal<-min(signallimits)
        maxsignal<-max(signallimits)
      }
      if(limsignals==T & levelz%in%comparisons){
        minsignal<-(-1) #theses are -1 and 1 per default. (but one could change them if necessary)
        maxsignal<-(1)
        if(class(signalcomplimits)!="NULL"){
          minsignal<-min(signalcomplimits)
          maxsignal<-max(signalcomplimits)
        }
      }
    }

    #final dataset per conditionlevel:
    datmat2<-interpdataset[[levelz]]
    binwidthlevel = 1
    if(levelz%in%c("comp1", "comp2", "comp3",  "comp4",  "comp5", "difference")){
      binwidthlevel = 0.1
    }
    # plot it -----
    whatsthat<-pretty(seq(minsignal, maxsignal, length.out =10), n = 10)
    labels2keep<-c(which.min(whatsthat), which(whatsthat==0), which.max(whatsthat))
    labelswhats<-whatsthat
    labelswhats[-labels2keep]<-""



    plot<-ggplot2::ggplot(datmat2, ggplot2::aes(x, y, z = value)) +
      ggplot2::geom_tile(ggplot2::aes(fill = value)) +
      ggplot2::geom_point(ggplot2::aes(x=x, y=y, fill=value, col=value), shape=15) + #, binwidth = 0.1
      {if (quick)ggplot2::geom_contour(colour = 'white', alpha = 0.5, size=0.3, breaks = whatsthat)}+
      {if (!quick)ggplot2::geom_contour(colour = 'white', alpha = 0.5, size=0.3, breaks = whatsthat)}+

      {if (quick)ggplot2::geom_path(data = circledat2, ggplot2::aes(x, y, z = NULL), size=1.5, col="white") }+
      {if (quick)ggplot2::geom_path(data = circledat3, ggplot2::aes(x, y, z = NULL), size=1.5, col="white") }+
      {if (!quick)ggplot2::geom_path(data = circledat2, ggplot2::aes(x, y, z = NULL), size=2, col="white") }+

      ggplot2::geom_path(data = circledat, ggplot2::aes(x, y, z = NULL), size=0.5) +
      # draw the nose (haven't drawn ears yet)
      ggplot2::geom_line(data = data.frame(x = c(0.45, 0.5, 0.55), y = c(1.0, 1.05, 1.0)),
                         ggplot2::aes(x, y, z = NULL), linewidth=0.6) +
      # add points for the electrodes
      {if (quick)ggplot2::geom_point(data = dat, ggplot2::aes(x, y, z = NULL, fill = NULL),
                            shape = 21, colour = 'black', size = 0.5) }+
      {if (!quick) ggplot2::geom_point(data = dat, ggplot2::aes(x, y, z = NULL, fill = NULL),
                              shape = 21, colour = 'black', size = 0.5)}+
      ggplot2::labs(fill = paste0(levelz), col = paste0(levelz)) +
      ggplot2::theme_bw()+
      ggplot2::coord_fixed()+
      ggplot2::theme(legend.position = "bottom",legend.direction = "horizontal",
            legend.key.size = ggplot2::unit(1, "cm"),
            text = ggplot2::element_text(size=18),
            axis.line=ggplot2::element_blank(),axis.text.x=ggplot2::element_blank(),
            axis.text.y=ggplot2::element_blank(),axis.ticks=ggplot2::element_blank(),
            axis.title.x=ggplot2::element_blank(),
            axis.title.y=ggplot2::element_blank(),
            panel.background=ggplot2::element_blank(),panel.border=ggplot2::element_blank(),panel.grid.major=ggplot2::element_blank(),
            panel.grid.minor=ggplot2::element_blank(),plot.background=ggplot2::element_blank(),
            plot.margin = ggplot2::margin(1,1,1,1))+
      ggplot2::ggtitle(plottitle)

    if(palette!="custom"){

      plot <- plot +
        ggplot2::scale_fill_distiller(palette = palette, na.value = NA,  limits=c(min(whatsthat),max(whatsthat)),
                             oob=scales::squish, breaks = whatsthat,labels = labelswhats)+
        ggplot2::scale_colour_distiller(palette = palette, na.value = NA, limits=c(min(whatsthat),max(whatsthat)),
                               oob=scales::squish, breaks = whatsthat,labels = labelswhats)
    }


    if(palette=="custom"){
      plot<-plot +
        ggplot2::scale_fill_gradientn(
          colors=c( "#053061", "#2166ac", "#4393c3", "#92c5de", "#d1e5f0", "#f7f7f7", "#fddbc7", "#f4a582", "#d6604d", "#b2182b", "#67001f"),
          values=scales::rescale(c(minsignal,0, maxsignal)),oob=scales::squish,
          limits=c(min(whatsthat),max(whatsthat)),
          guide = ggplot2::guide_colorbar(ticks.colour = "white"), breaks = whatsthat, labels = labelswhats)+
        ggplot2::scale_colour_gradientn(
          colors=c(  "#053061", "#2166ac", "#4393c3", "#92c5de", "#d1e5f0", "#f7f7f7", "#fddbc7", "#f4a582", "#d6604d", "#b2182b", "#67001f"),
          values=scales::rescale(c(minsignal,0, maxsignal)),oob=scales::squish,
          limits=c(min(whatsthat),max(whatsthat)),breaks = whatsthat, labels = labelswhats)
    }
    #shall the electrode labels be plotted?
    if(addlabels==T){
      plot<-plot+ggplot2::geom_text(inherit.aes = FALSE ,
                           data=electrodelocs,ggplot2::aes(x=x, y=y, label=channel), size=4, hjust=1, vjust=1 )
    }
    if(class(markchans)!="NULL"){
      chanstohighlight = subset(electrodelocs, channel %in% markchans)
      if (quick){
        plot<-plot+ggplot2::annotate("point",x=chanstohighlight$x, y=chanstohighlight$y, size=marksize, shape = 21, fill = chanmarkcol)
      }else{
        plot<-plot+ggplot2::annotate("point",x=chanstohighlight$x, y=chanstohighlight$y, size=marksize, shape = 21, fill = chanmarkcol)
      }

    }
    #set the legend accordingly for the comparisons (which ones)

    if(levelz%in%c("comp1", "comp2", "comp3", "comp4", "comp5")){

      if(levelz=="comp1"){
        plot<-plot+ ggplot2::labs(fill = paste0(comp1), col = paste0(comp1))
      }
      if(levelz=="comp2"){
        plot<-plot+ ggplot2::labs(fill = paste0(comp2), col = paste0(comp2))
      }
      if(levelz=="comp3"){
        plot<-plot+ ggplot2::labs(fill = paste0(comp3), col = paste0(comp3))
      }
      if(levelz=="comp4"){
        plot<-plot+ ggplot2::labs(fill = paste0(comp4), col = paste0(comp4))
      }
      if(levelz=="comp5"){
        plot<-plot+ ggplot2::labs(fill = paste0(comp5), col = paste0(comp5))
      }
    }

    #if only one difference plot (for 2 conditions)
    if(levelz%in%c("difference")){
      difference = paste0(abbrvlabel(conditionlevels[2])," - ",abbrvlabel(conditionlevels[1]))
      if(levelz=="difference"){
        plot<-plot+ ggplot2::labs(fill = paste0(difference), col = paste0(difference))
      }
    }
    plot<-plot+
      {if(!class(textfamily)=="NULL") ggplot2::theme(axis.title.x = ggplot2::element_text(family = textfamily))}+
      ggplot2::guides(fill = ggplot2::guide_colorbar(title.position = "top", title.hjust = 0.5))
    # remove legend
    if(extralegend){
      legendtitle = plot$labels$colour
      tempplot<-plot+
        {if(!class(textfamily)=="NULL") ggplot2::theme(text = ggplot2::element_text(family = textfamily))}+
        ggplot2::theme(legend.margin = ggplot2::margin(0,0,0,0, unit = "pt"),
              legend.key.size = ggplot2::unit(0.3, "cm"),
              legend.box.margin = ggplot2::margin(0,0,0,0, unit = "pt"),
              legend.title = ggplot2::element_text(size = textsize-2),
              legend.text = ggplot2::element_text(size = textsize-2,
                                         margin = ggplot2::margin(0,10,0,10, unit = "pt"))
        )
      leg <- ggpubr::get_legend(tempplot)
      leg<-ggpubr::as_ggplot(leg)
      #ggsave("testleg.png",leg, width=5,height=5, units ="cm")
      legendlist[[paste0(levelz, "_legend")]]<-leg

      plot<-plot+ggplot2::theme(legend.position = "none")+ggplot2::xlab(paste0(legendtitle))+
        ggplot2::theme(axis.title.x = ggplot2::element_text(size = textsize))

    }

    plotlist[[levelz]]<-plot #for each level add plot to list

  }

  #if there is only one condition level to plot just return the plot
  if(length(conditionlevels)!=1){
    if(length(legendlist)!=0){
      plotlist<-append(plotlist, legendlist)
    }
    return(plotlist)
  }else{
    return(plot)
  }
}



