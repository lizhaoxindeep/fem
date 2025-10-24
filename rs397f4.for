$debug

*  STRESS ANALYSIS WITH FINITE ELEMENT METHOD FOR SOLVING
*  ELASTO PLASTIC MODEL PROBLEM -ADVANCING FACE-

* nodes                          < 9998
* nodes *2                         19996
* elements                       < 9999
* band width for stress analysis < 998
* band width for other analysis  < 199
* nodes for boundary condition   < 999
*
* ept
*  100: ’èí“`”M‰ðÍ
*  102: ’èí”M‰ž—Í‰ðÍi•ÏˆÊS‘©ðŒ‚Ì‚ÝŠm”Fj
*  103Fü–c’£—¦‚Ìˆá‚¢‚É‚æ‚é”M‰ž—ÍŠm”FAŠî€‰·“xŠm”FA
*       ‰×d‹«ŠEðŒŠm”F
*  104: ”ñ’èí‰ði–¢Š®j
*  105F”ñ’èí‰ði‚Å‚«‚½j
*  106F”ñ’èí”M‰ž—Í‰ðÍi‚½‚Ô‚ñ‚Å‚«‚½A’e‘Y«‚Í–¢Šm”Fj
*  107: ”ñ’èí“`”M‰ðÍ‚Å‰·“xðŒ‚ªŽžŠÔ‚ÌŠÖ”‚É‚È‚éê‡

*  201: “€Œ‹–c’£‚Ìl—¶i’e‘Y«‚ª‚¨‚©‚µ‚­‚È‚Á‚Ä‚µ‚Ü‚Á‚½j
*  204: “€Œ‹–c’£‚Æ’e‘Y«‚Ì‹¤‘¶
*  205: “€Œ‹–c’£“™‚É‚æ‚é”M‰ž—Í‚É‹Nˆö‚·‚é’e‘Y«
*  206: ”M“`’B‚Ìl—¶
*  207: ”M•úŽË
*  208: —Z‰ð”M¥“úŽË
*  209: ’e‘Y«‚É‚·‚é‚©‚µ‚È‚¢‚©‚ÌƒXƒCƒbƒ`‚ð‚Â‚¯‚½
*       ’e«‚Ìê‡”j‰ó‚µ‚È‚¢‚½‚ß‚É•K—v‚ÈˆêŽ²ˆ³k‹­“x
*       ‚Æ”j‰ó‚µ‚È‚¢‚½‚ß‚É•K—v‚Èˆø’£‹­“x‚ðo—Í‚·‚é‚æ‚¤‚É‚µ‚½
*
* eptf
*  100: ’èí–O˜aZ“§—¬
*  101: ”ñ’èí–O˜aZ“§—¬
*  102: —¬‘¬‚ðŒvŽZ
*
* rs
*  107: ‚·‚×‚Äƒoƒ“ƒhƒ}ƒgƒŠƒNƒX‰»—v‘f”9999ŒÂB
*  111: •X‚Æ…‚Ì”ä”M“±“ü

*—Z‰ð”M the heat of fusion 79.7 cal/g
*“€Œ‹Žž‚Ìƒ„ƒ“ƒO—¦•Ï‰»
*“€Œ‹Žž‚Ì“§…ŒW”•Ï‰»
*ü–c’£—¦‚Ì”ñüŒ`«
*4:Z“§—¬‰·“x‘Î—¬
*5:“§…ŒW”‰·“xE‰ž—ÍˆË‘¶
*6FŠÔŒ„…“€Œ‹
*”M‘Î—¬
*’èí—¬
*



      COMMON/C1/KKK(8),NOD(9999,3),X(9998,2),LOAD(19996),
     >          JP(999),HENI(999)
      COMMON/C2/E(9999),V(9999),THICK(9999),EARTH(9999),VOLU(9999),
     >          FOME,
     1          PHI(9999),qu(9999),RPHI(9999),Rqu(9999),TNS(9999),         
     2          RTNS(9999),DIL(9999),STRF(3)                              
      COMMON/C3/ST(19996,998),XMK(6,6),XX(19996),BMAT(3,6,9999),              
     1          STRESS(3,6,9999),XXT(19996)                                
      COMMON/C4/ISTAT(9999),SIGM(3,9999),STRM(3,9999),IIII,LK,FOM,KM,      
     1          ORX(9999),ORY(9999),NA,NAL                                
      COMMON/C5/QE(2,9999),QS(6,9999),PRST(3,9999),NC	  
	  COMMON/C6/a0(9999),c0(9999),es0(9999),akk(9999),bkk(9999),
     >          dkk(9999),tkk(9999),az0(9999),fk01(9999),alf
      dimension ifrac(9999),tcond(9999),
     >          fk1(9998,199),ff1(9998),jpt(999),pt(999,9980),
     >          fnodetemp(9998),expancoeff(9999),
     >          eletemp(9999),elethermstrain(9999),
     >          c1(3,3,9999),dens(9999),fc1(9998,199),capa(9999),
     >          cx(9999),cy(9999),
     >          fe2(9999),
     >          jtt(99,2),fluidtemp(99,9980),flen3(99),
     >          aarea(9999),poro(9999),sun(99,9980),fk(9999),
     >          fkr(9999),jpwp(149),pwp(149,9980),ss(9999),
     >          fnodepres(9998),elepresr(9999),flowk(9998,199),
     >          flowc(9998,199),flowv(9999,4),strmr(3,9999),
     >          prstr1(9999),prstr2(9999),prstr3(9999),
     >          ufw2(9999),sr2(9999),eletempr(9999)
     >          ,fmyu(9999),fdens(9999),
     >          pe1(9999),ft1(9999),fb1(9999),
     >          pe2(9999),ft2(9999),fb2(9999),
     >          elepres(9999),stsv(9999),
     >          alpha(9999),aaalp(9999),bbalp(9999),ppstr(9999),
     >          tss(9999),tssr(9999),
     >          tload(19996),pload(19996),rload(19996),alphar(9999)

*      DOUBLE PRECISION ST
      REAL LOAD
      character(20) datfile,outfile,pltfile,nodfile

      open(9,file='kes.tmp')
      write(*,*)'a,c,es,ak,bk,dk,t,k0,ec,ncon'
      write(*,*)'Data file name?'
      read(*,*)datfile

      lc=len_trim(datfile)
      outfile=datfile(1:lc-3)
      outfile(lc-3:lc)='.out'
      pltfile=datfile(1:lc-3)
      pltfile(lc-3:lc)='.plt'
      nodfile=datfile(1:lc-3)
      nodfile(lc-3:lc)='.nod'

      open(5,file=datfile)
      write(*,*)'Opening  ',datfile
      open(6,file=outfile)
      write(*,*)'Opening  ',outfile
      open(7,file=pltfile)
      write(*,*)'Opening  ',pltfile
      open(8,file=nodfile)
      write(*,*)'Opening  ',nodfile

C *** READING INPUT DATA ***                                            
C                                                                       
      NA=0

      write(*,*)'reading data'

    4 continue
      CALL INPUTD(ifrac,nelem,dump,tcond,npoin,jpt,pt,npresst,
     >            expancoeff,index,reftemp,ntemp,dt,dens,capa,cx,cy,
     >            fe2,ntrans,transcoef,emissivity,jtt,
     >            fluidtemp,flen3,poro,fusionheat,sun,
     >            istr,itempe,iflow, 
     >            pe1,ft1,fb1,pe2,ft2,fb2,
     >            jpwp,pwp,ss,npwp,
     >            fpini,accy,jbw2,ufw2,sr2,neleex,ipxc1,ipxc2,
     >            cthre,cscale,alpha,aaalp,bbalp,ncon)

      write(*,*)'Initializing'

* Initialization
        if(na.eq.0)then
        do 2090 i=1,nelem
        eletemp(i)=reftemp
        elethermstrain(i)=0.
        elepresr(i)=0.
        elepres(i)=fpini
        prstr1(i)=0.
        prstr2(i)=0.
        prstr3(i)=0.
        ifrac(i)=0.
        alpha(i)=0.
          do 2093 j=1,4
          flowv(i,j)=0.
 2093     continue
          do 2092 j=1,3
          sigm(j,i)=0.
          strm(j,i)=0.
          strmr(j,i)=0.
          prst(j,i)=0.
 2092   continue
 2090   continue

        do 2091 i=1,npoin
        fnodetemp(i)=0.
        fnodepres(i)=0.
 2091   continue

        endif

*        if(istr.ne.0)then
*        endif

* Loop for transient solver
      do 7222 itemp=1,ntemp

        if(itemp.ne.1)then
        do 3198 i=1,npoin*2
        load(i)=rload(i)
 3198   continue
        endif

      do 3199 i=1,nelem
      elepresr(i)=elepres(i)
 3199 continue

* Updating alpha

      do 8929 i=1,nelem
      alphar(i)=alpha(i)
 8929 continue

        if(itemp.ne.1)then
        do 5222 i=1,nelem
        aaaalp=aaalp(i)+bbalp(i)*(0.25*prst(2,i)+0.75*prst(1,i))/1.e6
          if(aaaalp.gt.1)then
          alpha(i)=1.
          endif
          if(aaaalp.lt.0)then
          alpha(i)=0.
          endif
          if(aaaalp.ge.0.and.aaaalp.le.1)then
          alpha(i)=aaaalp
          endif
 5222   continue
        endif

      write(6,*)
      write(6,8216)itemp
 8216 format(1h ,'itemp=',i5)

***************************************
* Iteration to converge solution

*      do 6826 iicc=1,10

      do 2865 i=1,nelem

      et=eletemp(i)

      a1=alog10(pe1(i))
      a2=alog10(pe2(i))
      a=a1+(a1-a2)*(et-ft1(i))/(ft1(i)-ft2(i))
      b=fb1(i)+(fb1(i)-fb2(i))*(et-ft1(i))/(ft1(i)-ft2(i))
      esiga=(1.+v(i))*(sigm(1,i)+sigm(2,i))/3.+alpha(i)*elepres(i)
      fk(i)=10**(a-b*esiga/1.e6)
*      fk(i)=10**(a+b*esiga/1.e6)
      fkr(i)=fk(i)


        if(et.gt.0.)then
        fmyu(i) =1.782e-3-0.557e-4*et+1.005e-6*et**2
     >          -0.936e-8*et**3+0.338e-10*et**4

*        write(*,*)fmyu(i)

        fdens(i)=1000.+5.27e-2*et-7.56e-3*et**2
     >          +4.23e-5*et**3-1.344e-7*et**4
        fk(i)=fkr(i)



        else
        fmyu(i)=1.782e-3
        fdens(i)=1.000
          if(et.gt.fe2(i))then
          fk(i)=fkr(i)*(1.-et/fe2(i))
            if(fk(i).lt.fkr(i)/1000.)fk(i)=fkr(i)/1000.
          else
          fk(i)=fkr(i)/1000.
          endif
        endif
 2865 continue
      write(*,*)'itemp',itemp

C *** CALCULATION OF STIFFNESS MATRIX ***                               
C

*      write(*,*)itemp,na

        if(na.eq.0)then

* istr   0: No stress analysis
*	 1: Elastic stress analysis
*	 2: Elasto-plastic stress analysis

* itempe 0: No heat solver
*	 1: Steady state heat analysis
*	 2: Transient heat analysis
*        3: Transient heat analysis (variable boundary conditions)

* iflow  0: No fluid flow
*        1: Steady fluid flow
*        2: Transient fluid flow
*        3: Transient fluid flow (variable boundary conditions)

* Calculation on temperature

          if(itempe.ne.0)then

*********************************************************
* Iteration for freezing and thawing
      do 2444 ncyc=1,10

            call tstiff(nelem,npoin,fk1,tcond,fc1,
     >                  dens,capa,
     >                  aarea,jbw2,eletemp,fe2,poro,ufw2,sr2,
     >                  eletempr,fusionheat)
              if(itempe.eq.1)then
              call tanalys
     >        (fk1,ff1,npoin,jpt,pt,npresst,fnodetemp,jbw2)
              else
              call transient(fk1,fc1,npoin,fnodetemp,itemp,reftemp,
     >                       jpt,pt,dt,npresst,
     >                       ntrans,transcoef,fluidtemp,flen3,jtt,
     >                       emissivity,eletemp,
     >                       nelem,sun,jbw2,
     >                       istr,fnodepres,ncyc)
              endif

 2444 continue
*********************************************************

            call calavetemp(fnodetemp,nelem,eletemp)
 

              if((istr.eq.0.and.itempe.eq.2).or.
     >           (istr.eq.0.and.itempe.eq.3).or.
     >           (istr.eq.0.and.itempe.eq.1))then
                if(iflow.eq.0)then
*                write(*,*)'elemout for temp'
*                call elemout(cx,cy,strmr,prstr1,prstr2,prstr3,
*     >                       ifrac,flowv,nelem)
*		  if(ncyc.eq.10)then

*                  do 7287 i=1,npoin
*                  write(8,3237)I,fnodetemp(i),fnodepres(i)/1.e6
* 7287             continue
*                  endif
                goto7222
                endif
              endif

              if(istr.eq.0.and.itempe.eq.1)goto3


              if((istr.ge.1.and.istr.le.2).and.
     >           (itempe.ge.1.and.itempe.le.3))then
              call thermalforce
     >        (index,eletemp,nelem,expancoeff,reftemp,
     >         fe2,elethermstrain,ufw2,sr2,poro,tss,tssr,
     >         npoin,tload)
              else
              do 7259 i=1,npoin*2
              tload(i)=0.
 7259         continue
              endif

          endif

        endif
*------------------------------------------------------------
* Calculation on pore fluid flow

      do 8865 i=1,nelem

      et=eletemp(i)

      a1=alog10(pe1(i))
      a2=alog10(pe2(i))
      a=a1+(a1-a2)*(et-ft1(i))/(ft1(i)-ft2(i))
      b=fb1(i)+(fb1(i)-fb2(i))*(et-ft1(i))/(ft1(i)-ft2(i))
      esiga=(1.+v(i))*(sigm(1,i)+sigm(2,i))/3.+alpha(i)*elepres(i)
      fk(i)=10**(a-b*esiga/1.e6)
      fkr(i)=fk(i)


        if(et.gt.0.)then
        fmyu(i) =1.782e-3-0.557e-4*et+1.005e-6*et**2
     >          -0.936e-8*et**3+0.338e-10*et**4

*        write(*,*)fmyu(i)

        fdens(i)=1000.+5.27e-2*et-7.56e-3*et**2
     >          +4.23e-5*et**3-1.344e-7*et**4
        fk(i)=fkr(i)

        else
        fmyu(i)=1.782e-3
        fdens(i)=1.000
          if(et.gt.fe2(i))then
          fk(i)=fkr(i)*(1.-et/fe2(i))
            if(fk(i).lt.fkr(i)/1000.)fk(i)=fkr(i)/1000.
          else
          fk(i)=fkr(i)/1000.
          endif
        endif
 8865 continue

          if(iflow.eq.1)then

*          write(*,*)'iflow=1'

          call fstiff(nelem,npoin,flowk,fk,fmyu,flowc,
     >                fdens,ss,aarea,jbw2)
          call fanalys(flowk,bmat,npoin,jpwp,pwp,npwp,
     >                 fnodepres,accy,jbw2,fnodetemp)
          call flowspeed(nelem,flowv,fnodepres,fk,fmyu,accy,
     >                   fnodetemp)

*      do 2111 i=1,nelem
*      write(*,*)'flowv',i,(flowv(i,j),j=1,4)
* 2111 continue

          endif

          if(iflow.eq.2.or.iflow.eq.3)then
      write(*,*)'calling fstiff'
          call fstiff(nelem,npoin,flowk,fk,fmyu,flowc,
     >                fdens,ss,aarea,jbw2)
      write(*,*)'calling ftransient'
          call ftransient(flowk,flowc,npoin,fnodepres,itemp,fpini,
     >                    jpwp,pwp,dt,npwp,accy,
     >                    jbw2,fnodetemp)
      write(*,*)'calling flowspeed'
          call flowspeed(nelem,flowv,fnodepres,fk,fmyu,accy,
     >                   fnodetemp)
          endif

*      write(*,*)'istr,itempe',istr,itempe

*          if(istr.eq.0.and.itempe.eq.0)then
*          write(*,*)'elemout for flow'

*******************************
* Average pore pressure

      do 1225 i=1,nelem
      sum=0.
      do 1226 j=1,3
      sum=sum+fnodepres(nod(i,j))
 1226 continue
      elepres(i)=sum/3.
 1225 continue
******************************

        if(iflow.ne.0)then

        write(*,*)'ppstrain is called'

        call ppstrain
     >  (index,elepres,nelem,elepresr,alpha,alphar,ppstr,pload)
        else
        do 3719 i=1,npoin*2
        pload(i)=0.
 3719   continue
        endif
*------------------------------------------------------------
*      write(*,*)'before stifns OK'

* Structural analysis

       do 3901 icon=1,ncon

       write(*,*)'icon',icon

    1   CALL STIFNS(c1,itemp,icon,prstr1,prstr2,prstr3)

C          IN THIS SUBROUTINE CALL FEM(XE)                              

C *** SOLVING STIFFNESS EQUATIONS AND PRINTING RESULTS ***              

      do 9739 i=1,npoin*2
      rload(i)=load(i)
 9739 continue

    2   continue
        CALL ANALYS(ifrac,dump,elethermstrain,c1,npoin,istr,
     >              neleex,ipxc1,ipxc2,cthre,cscale,stsv,
     >              elepres,alpha,
     >              pload,tload,prstr1,prstr2,prstr3)

 3901 continue

* 6826 continue
*************************************************************

C         IN THIS SUBROUTINE CALL SOLVER                                
C                            CALL ANAROK                                
C                            CALL ANACOL                                

*        endif

*      WRITE(6,3488)
*      WRITE(6,3498)(i,ORX(I),ORY(I),i=1,nelem)
* 3488 FORMAT(1h ,'CENTER COORDINATE OF EACH TRIANGULAR ELEMENT')
* 3498 FORMAT(4(I5,2F7.1))

*      write(6,*)
*      WRITE(6,3487)                                                     
* 3487 FORMAT(1h ,'NODAL DISPLACEMENT')
*      do 1999 i=1,npoin                        
*      WRITE(6,3497)I,(XXT(2*(I-1)+J),J=1,2)
* 3497 format(1h ,i5,2g12.4)
* 1999 continue

*      write(6,*)'itemp, na',itemp,na

         if(itempe.eq.0.or.itempe.eq.1.or.istr.eq.0)then

         write(*,*)' na= ',na

           if(na.lt.0)goto3
*           if(na.eq.0.and.iflow.eq.0)goto3
           if(na.gt.0)goto4
         endif

          call elemout(cx,cy,strmr,prstr1,prstr2,prstr3,
     >                 ifrac,flowv,nelem,elepres,stsv,alpha)

      do 7288 i=1,npoin
      WRITE(6,3497)I,(XXT(2*(I-1)+J),J=1,2)
*      write(8,3237)I,(XXT(2*(I-1)+J)/float(iiii+1),J=1,2),
*     >             fnodetemp(i),fnodepres(i)/1.e6
      write(8,3237)I,(XXT(2*(I-1)+J),J=1,2),fnodetemp(i),fnodepres(i)
 7288 continue
 3237 FORMAT(I5,5g12.4)
 3497 format(1h ,i5,2g12.4)


*          do 7288 i=1,npoin
*          write(8,3237)I,fnodetemp(i),fnodepres(i)/1.e6
* 3237     format(i5,"  0  0",2g12.5)
* 7288     continue
*          write(*,*)'before goto7222 ok',itemp,na
*          goto7222
*          endif

 7222 continue

    3 continue

      close(9)

      STOP
      END
****************************************************
      subroutine solver2(aaa,bbb,xxx,npoin,jbw2,switch)
      dimension aaa(9998,199),bbb(9998),xxx(9998),xxxt(9998)
*      real*8 temp, sum

      write(6,*)'solver2'
      do 3244 i=1,npoin
      do 3244 j=1,jbw2
      write(6,3983)i,j,aaa(i,j)
 3983 format(1h ,2i5,g12.4)
 3244 continue

      do 3243 i=1,npoin
      write(6,3982)i,bbb(i)
 3982 format(1h ,i5,g12.4)
 3243 continue


        if(aaa(1,1).lt.0.)then
        do 3200 i=1,npoin
        do 3199 j=1,jbw2
        aaa(i,j)=aaa(i,j)*switch
 3199   continue
        bbb(i)=bbb(i)*switch
 3200   continue
        endif

*      do 3199 i=1,npoin
*      write(6,3982)i,bbb(i)
* 3982 format(1h ,i5,g12.4)
* 3199 continue

*‘½•ª‚±‚±‚Ü‚Å‡‚Á‚Ä‚éB

*      write(*,*)'jbw2,npoin',jbw2,npoin

      do 120 i=1,npoin
      ip=npoin-i+1
      if(ip.gt.jbw2) ip=jbw2
      do 120 j=1,ip
      jq=jbw2-j
      if(jq.gt.i-1) jq=i-1
*      write(*,*)'i,j',i,j
      sum=aaa(i,j)
      if(jq.eq.0) go to 122
      do 124 k=1,jq
      ik=i-k
      jk=j+k
*      write(*,*)'ik,k+1,ik,jk',ik,k+1,ik,jk
  124 sum=sum-aaa(ik,k+1)*aaa(ik,jk)
  122 if(j.ne.1)goto126
*      write(*,*)'i,j,sum',i,j,sum
      if(sum.le.0.0) go to 5000
      temp=1./sqrt(sum)
      aaa(i,j)=temp
      go to 120
  126 aaa(i,j)=sum*temp
  120 continue

* solve bandmatrix

 1000 do 140 i=1,npoin
      j=i-jbw2+1
      if(i+1.le.jbw2) j=1
      sum=bbb(i)
      if(i.eq.1) go to 140
      do 142 k=j,i-1
      ia=i-k+1
  142 sum=sum-aaa(k,ia)*xxx(k)
  140 xxx(i)=sum*aaa(i,1)
      do 150 i=1,npoin
      ii=npoin-i+1
      j=ii+jbw2-1
      if(j.gt.npoin) j=npoin
      sum=xxx(ii)
      if(i.eq.1) go to 149
      ib=ii+1
      do 152 k=ib,j
      ia=k-ii+1
  152 sum=sum-aaa(ii,ia)*xxx(k)
  149 xxx(ii)=sum*aaa(ii,1)
  150 xxxt(ii)=xxxt(ii)+xxx(ii)
    6 continue
      return
 5000 write(*,5001) i,j
 5001 format(1h ,'!!!!!!!!!!!!!!!!!!!!ATENTION!!!!!!!!!!!!!!!!!!!!',/,
     >       1h ,'SOLUTION IS NOT OBTAINED!',/,
     >       1h ,'Sum is negative value in solver2. i, j =',2i5)

      stop
      end
**********************************************************************
      subroutine elemout(cx,cy,strmr,prstr1,prstr2,prstr3,
     >                   ifrac,flowv,nelem,elepres,stsv,alpha)
      COMMON/C4/ISTAT(9999),SIGM(3,9999),STRM(3,9999),IIII,LK,FOM,KM,      
     1          ORX(9999),ORY(9999),NA,NAL                                
      COMMON/C5/QE(2,9999),QS(6,9999),PRST(3,9999),NC
*	  COMMON/C6/a0,c0,es0,d0,s0,al0,az0,fk01,z
      dimension flowv(9999,4),
     >          ifrac(9999),
     >          cx(9999),cy(9999),strmr(3,9999),
     >          prstr1(9999),prstr2(9999),prstr3(9999),
     >          elepres(9999),stsv(9999),alpha(9999)

      WRITE(6,3398)
      do 100 lk=1,nelem
*      WRITE(6,3399)LK,cx(lk),cy(lk),
*     >           ((STRM(I,LK)+strmr(i,lk))*1.e6,I=1,3),
*     >             prstr1(lk)*1.e6,prstr2(lk)*1.e6,prstr3(lk),
*     1            (SIGM(I,LK)/1.e6,I=1,3),
*     >             PRST(1,LK)/1.e6,prst(2,lk)/1.e6,prst(3,lk),
*     2             Ifrac(LK),(flowv(lk,i),i=1,4),elepres(lk)/1.e6
  100 continue

      do 200 lk=1,nelem
      WRITE(7,3399)LK,cx(lk),cy(lk),
     >           ((STRM(I,LK)+strmr(i,lk))*1.e6,I=1,3),
     >             prstr1(lk)*1.e6,prstr2(lk)*1.e6,prstr3(lk),
     >            (SIGM(I,LK)/1.e6,I=1,3),
     >             PRST(1,LK)/1.e6,prst(2,lk)/1.e6,prst(3,lk),
     2             Ifrac(LK),(flowv(lk,i),i=1,4),elepres(lk)/1.e6,
     >             stsv(lk),alpha(lk)
  200 continue

 3399 FORMAT(1H ,I5,2f10.3,5g12.3,7f10.3,I6,3g12.3,f10.2,2g12.4,f10.2)
 3398 FORMAT(1H ,'ELEM ',
     >          'center-x  center-y  EPSX      EPSY      GAMXY     ',
     >                              'eps1      eps2      angle     ',
     >          'SIGX      SIGY      TAUXY     ',
     >          'SIG1      SIG2      ANGLE     ifrac ',
     >          'flowvx      flowvy      flowv      theta     ')
      return
      end
********************************************************************
      subroutine flowspeed(nelem,flowv,fnodepres,fk,fmyu,accy,
     >                     fnodetemp)

      COMMON/C1/KKK(8),NOD(9999,3),X(9998,2),LOAD(19996),
     >          JP(999),HENI(999)
      real load
      COMMON/C2/E(9999),V(9999),THICK(9999),EARTH(9999),VOLU(9999),
     >          FOME,     
     1          PHI(9999),qu(9999),RPHI(9999),Rqu(9999),TNS(9999),         
     2          RTNS(9999),DIL(9999),STRF(3)                              

      dimension flowv(9999,4),fnodepres(9998),fk(9999)
     >          ,fmyu(9999),fnodetemp(9998)

      pai=3.141593

*      write(9,*)fnodepres(10)

      do 100 i=1,nelem

      n1=nod(i,1)
      x1=x(n1,1)
      y1=x(n1,2)

      n2=nod(i,2)
      x2=x(n2,1)
      y2=x(n2,2)

      n3=nod(i,3)
      x3=x(n3,1)
      y3=x(n3,2)

      a1=x2*y3-x3*y2
      a2=x3*y1-x1*y3
      a3=x1*y2-x2*y1

      area=(a1+a2+a3)/2.

      b1=y2-y3
      b2=y3-y1
      b3=y1-y2

      c1=x3-x2
      c2=x1-x3
      c3=x2-x1

      con=fk(i)/fmyu(i)/(area*2.)*thick(i)

*      write(*,*)i,fk(i),fmyu(i),area,thick(i)
 
       et=fnodetemp(n1)
          if(et.lt.0.)et=0.
        fd1=1000.+5.27e-2*et-7.56e-3*et**2
     >          +4.23e-5*et**3-1.344e-7*et**4
       et=fnodetemp(n2)
          if(et.lt.0.)et=0.
        fd2=1000.+5.27e-2*et-7.56e-3*et**2
     >          +4.23e-5*et**3-1.344e-7*et**4
       et=fnodetemp(n3)
          if(et.lt.0.)et=0.
        fd3=1000.+5.27e-2*et-7.56e-3*et**2
     >          +4.23e-5*et**3-1.344e-7*et**4

      fp1=fnodepres(n1)
*        if(fp1.lt.0.)fp1=0.
      fp2=fnodepres(n2)
*        if(fp2.lt.0.)fp2=0.
      fp3=fnodepres(n3)
*        if(fp3.lt.0.)fp3=0.

      pres1=(fp1+x(n1,2)*accy*fd1)
      pres2=(fp2+x(n2,2)*accy*fd2)
      pres3=(fp3+x(n3,2)*accy*fd3)

*      pres1=fnodepres(n1)
*      pres2=fnodepres(n2)
*      pres3=fnodepres(n3)

      vx=-(b1*pres1+b2*pres2+b3*pres3)*con
      vy=-(c1*pres1+c2*pres2+c3*pres3)*con
      va=sqrt(vx**2+vy**2)

*        if(va.gt.0.1)
*     >  write(12,1111)i,vx,vy,va,b1,b2,b3,c1,c2,c3,con,
*     > fk(i),fmyu(i),area,thick(i),pres1,pres2,pres3
* 1111   format(1h ,i5,17g12.4)


      flowv(i,1)=vx
      flowv(i,2)=vy
      flowv(i,3)=va

        if(vx.ne.0.)then
        theta=atan(vy/vx)*180./pai
          if(vx.lt.0)then
          theta=theta+180.
          endif
        endif

        if(vx.eq.0.and.vy.ge.0.)then
        theta=90.
        endif

        if(vx.eq.0.and.vy.lt.0.)then
        theta=-90.
        endif

      flowv(i,4)=theta

  100 continue

*      do 2111 i=1,nelem
*      write(*,*)'flowv',i,(flowv(i,j),j=1,4)
* 2111 continue

      return
      end
*********************************************************************
      subroutine ftransient(flowk,flowc,npoin,fnodepres,itemp,fpini,
     >                      jpwp,pwp,dt,npwp,accy,
     >                      jbw2,fnodetemp)

      COMMON/C1/KKK(8),NOD(9999,3),X(9998,2),LOAD(19996),
     >          JP(999),HENI(999)
      dimension flowk(9998,199),flowc(9998,199),fnodepres(9998),
     >          ta(9998,199),tb(9998),jpwp(149),pwp(149,9980)
     >          ,fnodetemp(9998)
      real load

        if(itemp.eq.1)then
        do 300 i=1,npoin
        fnodepres(i)=fpini
  300   continue
        endif

      do 500 i=1,npwp
      fnodepres(jpwp(i))=pwp(i,itemp)
  500 continue

*      write(6,*)(fnodetemp(i),i=1,6)

      do 2100 i=1,npoin
      tb(i)=0.
 2100 continue

*      write(6,3232)(tb(i),i=1,6)
* 3232 format(1h ,'tb',6g12.4)

      do 100 i=1,npoin
      do 100 j=1,npoin
***************************
        if(j.ge.i)then
        ii=i
        jj=j-i+1
          if(jj.gt.jbw2)goto100
        else
        ii=j
        jj=i-j+1
          if(jj.gt.jbw2)goto100
        endif
***************************
        if(itemp.eq.1)then
        et=fnodetemp(j)
          if(et.lt.0.)et=0.
        fd=1000.+5.27e-2*et-7.56e-3*et**2
     >          +4.23e-5*et**3-1.344e-7*et**4
        tb(i)=tb(i)+(-0.5*flowk(ii,jj)+flowc(ii,jj)/dt)
     >       *(fnodepres(j)+x(j,2)*fd*accy)
        else
        tb(i)=tb(i)+(-0.5*flowk(ii,jj)+flowc(ii,jj)/dt)
     >       *fnodepres(j)
        endif
*      write(6,3333)i,j,tb(i),fk1(i,j),fc1(i,j),fnodetemp(j),dt
  100 continue

*      do 3200 i=1,npoin
*      do 3200 j=1,jbw2
*      write(6,3983)i,j,flowk(i,j),flowc(i,j)
* 3983 format(1h ,2i5,2g12.4)
* 3200 continue

      do 200 i=1,npoin
      do 200 j=1,jbw2
      ta(i,j)=0.5*flowk(i,j)+flowc(i,j)/dt
*      write(6,3984)i,j,ta(i,j)
* 3984 format(1h ,2i5,g12.4)
  200 continue

*      do 3200 i=1,npoin
*      do 3200 j=1,jbw2
*      write(6,3983)i,j,ta(i,j)
* 3983 format(1h ,2i5,g12.4)
* 3200 continue
*‚±‚±‚Ìta‚Í‡‚Á‚Ä‚Ü‚·B

      do 6500 i=1,npwp
      jjpt=jpwp(i)
      ppt=pwp(i,itemp)
      do 6600 j=1,npoin
      jjjpt=jjpt-j+1
*        if(jjjpt.le.0)goto6600
        if(jjjpt.ge.1)then
        jj=j
        jjj=jjjpt
        else
        jj=jjpt
        jjj=j-jjpt+1
        endif

*      write(11,*)'i,j,jj,jjj',i,j,jj,jjj

      if(jjj.le.jbw2)then
      tta=ta(jj,jjj)
      tb(j)=tb(j)-tta*ppt
      endif
 6600 continue
      ta(jjpt,1)=-1.e20
 6500 continue

*      do 3200 i=1,npoin
*      do 3200 j=1,jbw2
*      write(6,3983)i,j,ta(i,j)
* 3983 format(1h ,2i5,g12.4)
* 3200 continue
*‚±‚±‚Ìta‚à‚ ‚Á‚Ä‚Ü‚·Bsolver‚Ì’†‚©H
*      call sweep(ta,tb,fnodepres,1999,npoin)
      call solver2(ta,tb,fnodepres,npoin,jbw2,-1.)

      do 5500 i=1,npwp
      fnodepres(jpwp(i))=fnodepres(jpwp(i))+pwp(i,itemp)
 5500 continue

      write(6,*)' '
      write(6,*)'node, x, y, pore pressure (MPa)'
*      do 8999 i=1,npoin
*      write(6,9001)i,x(i,1),x(i,2),fnodepres(i)/1.e6
 9001 format(1h ,i5,3f10.2)
* 8999 continue

      return
      end
*************************************************************
      subroutine fstiff(nelem,npoin,flowk,fk,fmyu,flowc,
     >                  fdens,ss,aarea,jbw2)

      COMMON/C1/KKK(8),NOD(9999,3),X(9998,2),LOAD(19996),
     >          JP(999),HENI(999)
      real load
      dimension flowk(9998,199),fk(9999),ss(9999),
     >          flowc(9998,199),aarea(9999)
     >          ,fmyu(9999),fdens(9999)

      do 200 i=1,npoin
      do 200 j=1,jbw2
      flowk(i,j)=0.
      flowc(i,j)=0.
  200 continue

      do 100 i=1,nelem

      n1=nod(i,1)
      x1=x(n1,1)
      y1=x(n1,2)

      n2=nod(i,2)
      x2=x(n2,1)
      y2=x(n2,2)

      n3=nod(i,3)
      x3=x(n3,1)
      y3=x(n3,2)

      a1=x2*y3-x3*y2
      a2=x3*y1-x1*y3
      a3=x1*y2-x2*y1

      area=(a1+a2+a3)/2.
      aarea(i)=area

      b1=y2-y3
      b2=y3-y1
      b3=y1-y2

      c1=x3-x2
      c2=x1-x3
      c3=x2-x1

* Magnification to prevent error in SWEEP

      con=fk(i)/(fmyu(i)*4.*area)

      flowk(n1,1)=flowk(n1,1)+con*(b1*b1+c1*c1)
      if(n2.ge.n1)flowk(n1,n2-n1+1)=flowk(n1,n2-n1+1)+con*(b1*b2+c1*c2)
      if(n3.ge.n1)flowk(n1,n3-n1+1)=flowk(n1,n3-n1+1)+con*(b1*b3+c1*c3)

      if(n1.ge.n2)flowk(n2,n1-n2+1)=flowk(n2,n1-n2+1)+con*(b1*b2+c1*c2)
      flowk(n2,1)=flowk(n2,1)+con*(b2*b2+c2*c2)
      if(n3.ge.n2)flowk(n2,n3-n2+1)=flowk(n2,n3-n2+1)+con*(b2*b3+c2*c3)

      if(n1.ge.n3)flowk(n3,n1-n3+1)=flowk(n3,n1-n3+1)+con*(b1*b3+c1*c3)
      if(n2.ge.n3)flowk(n3,n2-n3+1)=flowk(n3,n2-n3+1)+con*(b2*b3+c2*c3)
      flowk(n3,1)=flowk(n3,1)+con*(b3*b3+c3*c3)

* Magnification to prevent error in SWEEP

      con2=ss(i)/(fdens(i)*9.81)*area/12.

*      call writeflowc(flowc)

      flowc(n1,1)=flowc(n1,1)+con2*2.
      if(n2.ge.n1)flowc(n1,n2-n1+1)=flowc(n1,n2-n1+1)+con2
      if(n3.ge.n1)flowc(n1,n3-n1+1)=flowc(n1,n3-n1+1)+con2
*      call writeflowc(flowc)


      if(n1.ge.n2)flowc(n2,n1-n2+1)=flowc(n2,n1-n2+1)+con2
      flowc(n2,1)=flowc(n2,1)+con2*2.
      if(n3.ge.n2)flowc(n2,n3-n2+1)=flowc(n2,n3-n2+1)+con2
*      call writeflowc(flowc)

      if(n1.ge.n3)flowc(n3,n1-n3+1)=flowc(n3,n1-n3+1)+con2
      if(n2.ge.n3)flowc(n3,n2-n3+1)=flowc(n3,n2-n3+1)+con2
      flowc(n3,1)=flowc(n3,1)+con2*2.
*      call writeflowc(flowc)

  100 continue

*Introduction of heat transfer
*
*        if(ntrans.ne.0)then
*        do 1200 i=1,ntrans
*        ip1=jtt(i,1)
*        ip2=jtt(i,2)
*        flowk(ip1,ip1)=flowk(ip1,ip1)+transcoef*flen3(i)/3.
*        flowk(ip2,ip2)=flowk(ip2,ip2)+transcoef*flen3(i)/3.
*        flowk(ip1,ip2)=flowk(ip1,ip2)+transcoef*flen3(i)/6.
*        flowk(ip2,ip1)=flowk(ip2,ip1)+transcoef*flen3(i)/6.
* 1200   continue
*        endif

*      do 6700 i=1,npoin
*      do 6700 j=1,jbw2
*      write(6,9990)i,j,flowc(i,j)
* 6700 continue
* 9990 format(2i5,g12.4)

      return
      end

*************************************************************************

*      subroutine writeflowc
*      dimension flowc(9998,1998
*      do 6800 i=1,4
*      write(6,9990)i,(flowc(i,j),j=1,4)
* 6800 continue
* 9990 format(i5,6g12.4)
*      return
*      end

      subroutine fanalys(flowk,bmat,npoin,jpwp,pwp,npwp,
     >                   fnodepres,accy,jbw2,fnodetemp)
      COMMON/C1/KKK(8),NOD(9999,3),X(9998,2),LOAD(19996),
     >          JP(999),HENI(999)
      real load
      dimension flowk(9998,199),bmat(9998),jpwp(149),pwp(149,9980),
     >          fnodepres(9998)
     >          ,fnodetemp(9998)

*      write(*,*)accy,fdens

      do 402 i=1,npoin
      bmat(i)=0.
  402 continue

      do 400 i=1,npoin
      do 400 j=1,npoin
        if(j.ge.i)then
        ii=i
        jj=j-i+1
          if(jj.gt.jbw2)goto400
        else
        ii=j
        jj=i-j+1
          if(jj.gt.jbw2)goto400
        endif

        et=fnodetemp(i)
          if(et.lt.0.)et=0.
        fd=1000.+5.27e-2*et-7.56e-3*et**2
     >       +4.23e-5*et**3-1.344e-7*et**4
      bmat(j)=bmat(j)-flowk(ii,jj)*x(i,2)*accy*fd
*????????????????????????????????????????????????????
  400 continue

      do 500 i=1,npwp
      jjpwp=jpwp(i)
      ppwp=pwp(i,1)
*      write(*,*)i,npwp,jjpwp,ppwp
      do 600 j=1,npoin
      jjjpwp=jjpwp-j+1
*        if(jjjpwp.le.0)goto600
        if(jjjpwp.ge.1)then
        jj=j
        jjj=jjjpwp
        else
        jj=jjpwp
        jjj=j-jjpwp+1
        endif

*      write(*,*)'jj,jjj',jj,jjj

      if(jjj.le.jbw2)then
      ffk1=flowk(jj,jjj)
      bmat(j)=bmat(j)-ffk1*ppwp
      endif
  600 continue
      flowk(jjpwp,1)=sign(1.e20,flowk(1,1))
  500 continue

*      do 700 i=1,6
*      write(6,9990)ff1(i),(fk1(i,j),j=1,6)
*  700 continue
* 9990 format(1h ,7g12.4)

*      CALL SWEEP(flowk,bmat,fnodepres,1999,Npoin)
      call solver2(flowk,bmat,fnodepres,npoin,jbw2,-1.)

      do 100 i=1,npwp
      jjpt=jpwp(i)
      fnodepres(jjpt)=fnodepres(jjpt)+pwp(i,1)
  100 continue

*      write(6,*)' '
*      write(6,*)'node, x, y, pore pressure (MPa)'
*      do 8999 i=1,npoin
*      write(6,9001)i,x(i,1),x(i,2),fnodepres(i)/1.e6
 9001 format(1h ,i5,3f10.2)
* 8999 continue
      return
      end

******************************************************************

      subroutine transient(fk1,fc1,npoin,fnodetemp,itemp,reftemp,jpt,
     >                     pt,dt,npresst,
     >                     ntrans,transcoef,fluidtemp,flen3,jtt,
     >                     emissivity,eletemp,
     >                     nelem,sun,jbw2,
     >                     istr,fnodepres,ncyc)

      COMMON/C1/KKK(8),NOD(9999,3),X(9998,2),LOAD(19996),
     >          JP(999),HENI(999)
      real load
      dimension fk1(9998,199),fc1(9998,199),fnodetemp(9998),
     >          ta(9998,199),tb(9998),jpt(999),pt(999,9980),
     >          fluidtemp(99,9980),flen3(99),jtt(99,2),
     >          eletemp(9999),eletempr(9999),
     >          sun(99,9980),fnodepres(9998),
     >          fnodetempr(9998)

        if(itemp.eq.1)then
        do 300 i=1,npoin
        fnodetemp(i)=reftemp
  300   continue
        endif

      do 500 i=1,npresst
      fnodetemp(jpt(i))=pt(i,itemp)
  500 continue

*ncyc‚ª1‚Ì‚Æ‚«‚Í‰Šú’l‚ð•Û‘¶A‚»‚¤‚Å‚È‚¢‚Æ‚«‚Í‰Šú’l‚É–ß‚·

        if(ncyc.eq.1)then
        do 611 i=1,npoin
        fnodetempr(i)=fnodetemp(i)
  611   continue
        else
        do 612 i=1,npoin
        fnodetemp(i)=fnodetempr(i)
  612   continue
        endif

*      write(6,*)(fnodetemp(i),i=1,6)

      do 2100 i=1,npoin
      tb(i)=0.
 2100 continue

*      write(6,3232)(tb(i),i=1,6)
* 3232 format(1h ,'tb',6g12.4)

      do 100 i=1,npoin
      do 100 j=1,npoin
***************************
        if(j.ge.i)then
        ii=i
        jj=j-i+1
          if(jj.gt.jbw2)goto100
        else
        ii=j
        jj=i-j+1
          if(jj.gt.jbw2)goto100
        endif
***************************
      tb(i)=tb(i)+(-0.5*fk1(ii,jj)+fc1(ii,jj)/dt)*fnodetemp(j)
*      write(6,3333)i,j,tb(i),fk1(i,j),fc1(i,j),fnodetemp(j),dt
  100 continue
 3333 format(1h ,2i5,5g12.4)

*Introduction of heat transfer and radiation

        if(ntrans.ne.0)then
        do 2001 i=1,ntrans
        ip1=jtt(i,1)
        ip2=jtt(i,2)
        avetemp=(fnodetemp(ip1)+fnodetemp(ip2))/2.+273.15

*      write(6,*)avetemp

        qrad=5.67*emissivity*(avetemp/100.)**4
*        if(itemp.eq.1.and.i.eq.1)write(10,*)'before',tb(ip1)
*        write(10,*)'before',itemp,tb(ip1)
        tb(ip1)=tb(ip1)+transcoef
     >         *(fluidtemp(i,itemp)-fnodetemp(ip1))*flen3(i)/2.
     >         -qrad*flen3(i)/2.+sun(i,itemp)*flen3(i)/2.
*        if(itemp.eq.1.and.i.eq.1)write(10,*)'after',tb(ip1)
*        write(10,*)'after',itemp,tb(ip1)
        tb(ip2)=tb(ip2)+transcoef
     >         *(fluidtemp(i,itemp)-fnodetemp(ip2))*flen3(i)/2.
     >         -qrad*flen3(i)/2.+sun(i,itemp)*flen3(i)/2.
 2001   continue
        endif


* Initializing and saving element temperature

       if(itemp.eq.1)then
       do 3002 i=1,nelem
       eletempr(i)=reftemp
 3002  continue
       endif

      do 200 i=1,npoin
      do 200 j=1,jbw2
      ta(i,j)=0.5*fk1(i,j)+fc1(i,j)/dt
  200 continue

*      write(6,*)'dt,npoin',dt,npoin

*      write(*,*)
*      write(6,*)'i,fk1'
*      do 6710 i=1,6
*      write(6,9991)i,(fk1(i,j),j=1,6)
* 6710 continue
* 9991 format(i5,6g12.4)

*      write(*,*)
*      write(6,*)'i,fc1'
*      do 6712 i=1,6
*      write(6,9992)i,(fc1(i,j),j=1,6)
* 6712 continue
* 9992 format(i5,6g12.4)

*      write(*,*)' '
*      write(6,*)'i,tb,ta before adjust'
*      do 6703 i=1,6
*      write(6,9993)i,tb(i),(ta(i,j),j=1,6)
* 6703 continue
* 9993 format(i5,7g12.4)


      do 6500 i=1,npresst
      jjpt=jpt(i)
      ppt=pt(i,itemp)
      do 6600 j=1,npoin
      jjjpt=jjpt-j+1
        if(jjjpt.ge.1)then
        jj=j
        jjj=jjjpt
        else
        jj=jjpt
        jjj=j-jjpt+1
        endif

*kes5‚ÅƒGƒ‰[Bƒoƒ“ƒh‚É‚µ‚½‚Æ‚±‚ëB

*      write(*,*)'i,j,jj,jjj,jjpt,jjjpt',i,j,jj,jjj,jjpt,jjjpt

        if(jjj.le.jbw2)then
        tta=ta(jj,jjj)
        tb(j)=tb(j)-tta*ppt
        endif
 6600 continue
      ta(jjpt,1)=1.e20
 6500 continue

*       write(6,*)'i,tb,ta'
*      do 6700 i=1,6
*      write(6,9990)i,tb(i),(ta(i,j),j=1,6)
* 6700 continue
* 9990 format(i5,7g12.4)

*      call sweep(ta,tb,fnodetemp,9998,npoin)
      call solver2(ta,tb,fnodetemp,npoin,jbw2,1.)

      do 5500 i=1,npresst
      fnodetemp(jpt(i))=fnodetemp(jpt(i))+pt(i,itemp)
 5500 continue

      do 2160 i=1,nelem
      sum=0.
      do 2260 j=1,3
      sum=sum+fnodetemp(nod(i,j))
 2260 continue
      eletemp(i)=sum/3.
 2160 continue

***************************************************

      do 9861 i=1,nelem
      sum=0.
      do 9860 j=1,3
      sum=sum+fnodetemp(nod(i,j))
 9860 continue
      eletemp(i)=sum/3.
 9861 continue

       do 3003 i=1,nelem
       eletempr(i)=eletemp(i)
 3003  continue

*      write(6,*)' '
*      write(6,*)'node, x, y, temperature'
*      do 8999 i=1,npoin
*      write(6,9001)i,x(i,1),x(i,2),fnodetemp(i)
 9001 format(1h ,i5,3f10.2)
* 8999 continue

        if(istr.eq.0)then
          if(ncyc.eq.10)then
          do 8998 i=1,npoin
          write(8,3237)i,fnodetemp(i),fnodepres(i)/1.e6
 8998     continue
          endif
 3237   FORMAT(I5," 0 0",2g12.4)
        endif

      return
      end





******************************************************************
      SUBROUTINE ppstrain
     >(index,elepres,nelem,elepresr,alpha,alphar,ppstr,pload)
      COMMON/C1/KKK(8),NOD(9999,3),X(9998,2),LOAD(19996),
     >          JP(999),HENI(999) 
      COMMON/C2/E(9999),V(9999),THICK(9999),EARTH(9999),VOLU(9999),
     >          FOME,     
     1          PHI(9999),qu(9999),RPHI(9999),Rqu(9999),TNS(9999),         
     2          RTNS(9999),DIL(9999),STRF(3)                              
      dimension elepres(9999),alpha(9999),ppstr(9999),elepresr(9999),
     >          pload(19996),alphar(9999)

      real load

*      do 400 i=1,npoin*2
*      load(i)=0.
*  400 continue

*
C strain by pp
*
      do 200 i=1,nelem

      bulk=e(i)/(3.*(1.-2.*v(i)))

* Linear strain by pp
          if(i.eq.1)then
          write(*,*)'elepres, elepresr',elepres(1),elepresr(1)
          endif
*        ppvolstr=alpha(i)*elepres(i)/bulk
        ppvolstr=(alpha(i)*elepres(i)-alphar(i)*elepresr(i))/bulk
        pplinstr=ppvolstr/3.
* Plane strain
        if(index.eq.0)then
        ppstr(i)=pplinstr*(1.+v(i))
        con=e(i)*thick(i)*pplinstr/(2.*(1.-2.*v(i)))
        else
* Plane stress
        ppstr(i)=pplinstr
        con=e(i)*thick(i)*pplinstr/(2.*(1.-v(i)))
        endif

      x1=x(nod(i,1),1)
      x2=x(nod(i,2),1)
      x3=x(nod(i,3),1)
      y1=x(nod(i,1),2)
      y2=x(nod(i,2),2)
      y3=x(nod(i,3),2)

      b1=y2-y3
      b2=y3-y1
      b3=y1-y2
      c1=x3-x2
      c2=x1-x3
      c3=x2-x1

      pload(nod(i,1)*2-1)=pload(nod(i,1)*2-1)+con*b1
      pload(nod(i,1)*2  )=pload(nod(i,1)*2  )+con*c1
      pload(nod(i,2)*2-1)=pload(nod(i,2)*2-1)+con*b2
      pload(nod(i,2)*2  )=pload(nod(i,2)*2  )+con*c2
      pload(nod(i,3)*2-1)=pload(nod(i,3)*2-1)+con*b3
      pload(nod(i,3)*2  )=pload(nod(i,3)*2  )+con*c3

  200 continue

*      write(*,*)' '
*      write(*,*)'Nodal force due to pp expansion'
*      do 300 i=1,npoin*2
*      write(*,1000)i,pload(i)
 1000 format(1h ,i5,g12.4)
*  300 continue

      return
      end
******************************************************************

      SUBROUTINE thermalforce
     >(index,eletemp,nelem,expancoeff,reftemp,
     > fe2,elethermstrain,ufw2,sr2,poro,tss,tssr,npoin,
     > tload)
      COMMON/C1/KKK(8),NOD(9999,3),X(9998,2),LOAD(19996),
     >          JP(999),HENI(999) 
      COMMON/C2/E(9999),V(9999),THICK(9999),EARTH(9999),VOLU(9999),
     >          FOME,     
     1          PHI(9999),qu(9999),RPHI(9999),Rqu(9999),TNS(9999),         
     2          RTNS(9999),DIL(9999),STRF(3)                              
      dimension expancoeff(9999),eletemp(9999),
     >          fe2(9999),
     >          elethermstrain(9999),ufw2(9999),sr2(9999),poro(9999),
     >          tss(9999),tssr(9999),tload(19996)

      real load

      do 400 i=1,npoin*2
      tload(i)=0.
  400 continue

*
C thermal force
*
      do 200 i=1,nelem


        if(reftemp.gt.fe2(i).and.reftemp.lt.0.)then
        write(*,*)'reftemp should not be between fe2 and less than 0.'
        return
        endif

*        if(index.eq.0)then
*        fsdevcoeff=2.*(1+v(i))
*        else
        fsdevcoeff=3.
*        endif

* Linear strain of water or ice at eletemp
        if(eletemp(i).gt.0.)then
        flsw=0.000137237-5.22526e-5*eletemp(i)
     >      +7.51338e-6*eletemp(i)**2
     >      -4.05718e-8*eletemp(i)**3+1.40139e-10*eletemp(i)**4
        else
        flsw=0.0906103+6.40993e-5*eletemp(i)+1.48253e-7*eletemp(i)**2
          if(eletemp(i).gt.fe2(i))then
          flsw=(flsw-0.000137237)*eletemp(i)/fe2(i)+0.000137237
          endif
        endif
      flsw=flsw*poro(i)/fsdevcoeff*(1.-ufw2(i))*sr2(i)
*Linear strain of water or ice at reftemp
        if(reftemp.gt.0.)then
        flswr=0.000137237-5.22526e-5*reftemp+7.51338e-6*reftemp**2
     >   -4.05718e-8*reftemp**3+1.40139e-10*reftemp**4
        else
        flswr=0.0906103+6.40993e-5*reftemp+1.48253e-7*reftemp**2
        endif
      flswr=flswr*poro(i)/fsdevcoeff*(1-ufw2(i))*sr2(i)
*Strain change due to water or ice
      dflsw=flsw-flswr
*
      tss(i)=(eletemp(i)-reftemp)*expancoeff(i)+dflsw
*        if(itemp.eq.1)then
        ts=tss(i)
*        else
*        ts=tss(i)-tssr(i)
*        endif
      tssr(i)=tss(i)

*      write(*,*)'ts',ts

* Plane strain
        if(index.eq.0)then
        elethermstrain(i)=ts*(1.+v(i))
        con=e(i)*thick(i)*ts/(2.*(1.-2.*v(i)))
        else
* Plane stress
        elethermstrain(i)=ts
        con=e(i)*thick(i)*ts/(2.*(1.-v(i)))
        endif

      x1=x(nod(i,1),1)
      x2=x(nod(i,2),1)
      x3=x(nod(i,3),1)
      y1=x(nod(i,1),2)
      y2=x(nod(i,2),2)
      y3=x(nod(i,3),2)

      b1=y2-y3
      b2=y3-y1
      b3=y1-y2
      c1=x3-x2
      c2=x1-x3
      c3=x2-x1

      tload(nod(i,1)*2-1)=tload(nod(i,1)*2-1)+con*b1
      tload(nod(i,1)*2  )=tload(nod(i,1)*2  )+con*c1
      tload(nod(i,2)*2-1)=tload(nod(i,2)*2-1)+con*b2
      tload(nod(i,2)*2  )=tload(nod(i,2)*2  )+con*c2
      tload(nod(i,3)*2-1)=tload(nod(i,3)*2-1)+con*b3
      tload(nod(i,3)*2  )=tload(nod(i,3)*2  )+con*c3

  200 continue

      write(6,*)' '
      write(6,*)'Nodal force due to thermal expansion'
      do 300 i=1,npoin*2
      write(6,1000)i,tload(i)
 1000 format(1h ,i5,g12.4)
  300 continue

      return
      end

***********************************************************************

      subroutine calavetemp
     >(fnodetemp,nelem,eletemp)
      COMMON/C1/KKK(8),NOD(9999,3),X(9998,2),LOAD(19996),
     >          JP(999),HENI(999)
      real load
      COMMON/C2/E(9999),V(9999),THICK(9999),EARTH(9999),VOLU(9999),
     >          FOME,
     1          PHI(9999),qu(9999),RPHI(9999),Rqu(9999),TNS(9999),         
     2          RTNS(9999),DIL(9999),STRF(3)                              
      dimension fnodetemp(9998),
     >          eletemp(9999)

*      write(6,*)' '
*      write(6,*)'Element, ox, oy, temperature'

      do 100 i=1,nelem
      sum=0.
*      sumx=0.
*      sumy=0.
      do 200 j=1,3
      sum=sum+fnodetemp(nod(i,j))
*      sumx=sumx+x(nod(i,j),1)
*      sumy=sumy+x(nod(i,j),2)
  200 continue
      eletemp(i)=sum/3.
*      write(6,1000)i,cx(i),cy(i),eletemp(i)
*        if(index.eq.0)then
*        elethermstrain(i)=(eletemp(i)-reftemp)*expancoeff(i)*(1.+v(i))
*        else
*        elethermstrain(i)=(eletemp(i)-reftemp)*expancoeff(i)
*        endif
  100 continue
 1000 format(1h ,i5,3f10.2)
      return
      end

***********************************************************************

      subroutine tstiff(nelem,npoin,fk1,tcond,fc1,
     >                  dens,capa,
     >                  aarea,jbw2,eletemp,fe2,poro,
     >                  ufw2,sr2,eletempr,fusionheat)

      COMMON/C1/KKK(8),NOD(9999,3),X(9998,2),LOAD(19996),
     >          JP(999),HENI(999)
      real load
      dimension fk1(9998,199),tcond(9999),capa(9999),
     >          fc1(9998,199),dens(9999),
     >          aarea(9999),eletemp(9999),
     >          fe2(9999),poro(9999),ufw2(9999),sr2(9999),
     >          eletempr(9999)

      do 200 i=1,npoin
      do 200 j=1,jbw2
      fk1(i,j)=0.
      fc1(i,j)=0.
  200 continue

      do 100 i=1,nelem

      n1=nod(i,1)
      x1=x(n1,1)
      y1=x(n1,2)

      n2=nod(i,2)
      x2=x(n2,1)
      y2=x(n2,2)

      n3=nod(i,3)
      x3=x(n3,1)
      y3=x(n3,2)

      a1=x2*y3-x3*y2
      a2=x3*y1-x1*y3
      a3=x1*y2-x2*y1

      area=(a1+a2+a3)/2.
      aarea(i)=area

      b1=y2-y3
      b2=y3-y1
      b3=y1-y2

      c1=x3-x2
      c2=x1-x3
      c3=x2-x1

********************

      avelt=(eletemp(i)+eletempr(i))/2.
      ratioiw1=avelt/fe2(i)
        if(ratioiw1.gt.1.)ratioiw1=1.
        if(ratioiw1.lt.0.)ratioiw1=0.

      ratioiw2=(1.-ufw2(i))*ratioiw1

      condi=2.2
      condw=0.5765+1.9852e-4*avelt+3.890e-5*avelt**2
     >     -4.271e-7*avelt**3+1.2686e-9*avelt**4
      tcondiw=tcond(i)
     >       +poro(i)*sr2(i)*condw*(1.-ratioiw2)
     >       +poro(i)*sr2(i)*condi*ratioiw2
**********************
      ratioiwa=eletemp(i)/fe2(i)
        if(ratioiwa.gt.1.)ratioiw1=1.
        if(ratioiwa.lt.0.)ratioiw1=0.
      ratioiwb=eletempr(i)/fe2(i)
        if(ratioiwb.gt.1.)ratioiw1=1.
        if(ratioiwb.lt.0.)ratioiw1=0.

      capai=2103.+3.186*avelt-0.03149*avelt**2
     >     -5.193e-5*avelt**3
      capaw=4200.
      capaiw=capa(i)
     >      +poro(i)*sr2(i)*capaw*(1.-ratioiw2)
     >      +poro(i)*sr2(i)*capai*ratioiw2

         if(eletemp(i)-eletempr(i).ne.0.)then
         capaiw=(capaiw*dens(i)*abs(eletemp(i)-eletempr(i))
     >         +fusionheat*poro(i)*sr2(i)*(abs(ratioiwa-ratioiwb)))
     >         /dens(i)/abs(eletemp(i)-eletempr(i))
         endif

       if(capaiw.le.0.)then
       capaiw=0.1
       write(*,*)'capaiw was set at 0.1'
       endif
********

      con=tcondiw/4./area

      fk1(n1,1)=fk1(n1,1)+con*(b1*b1+c1*c1)
      if(n2.ge.n1)fk1(n1,n2-n1+1)=fk1(n1,n2-n1+1)+con*(b1*b2+c1*c2)
      if(n3.ge.n1)fk1(n1,n3-n1+1)=fk1(n1,n3-n1+1)+con*(b1*b3+c1*c3)

      if(n1.ge.n2)fk1(n2,n1-n2+1)=fk1(n2,n1-n2+1)+con*(b1*b2+c1*c2)
      fk1(n2,1)=fk1(n2,1)+con*(b2*b2+c2*c2)
      if(n3.ge.n2)fk1(n2,n3-n2+1)=fk1(n2,n3-n2+1)+con*(b2*b3+c2*c3)

      if(n1.ge.n3)fk1(n3,n1-n3+1)=fk1(n3,n1-n3+1)+con*(b1*b3+c1*c3)
      if(n2.ge.n3)fk1(n3,n2-n3+1)=fk1(n3,n2-n3+1)+con*(b2*b3+c2*c3)
      fk1(n3,1)=fk1(n3,1)+con*(b3*b3+c3*c3)

      con2=dens(i)*capaiw*area/12.

      fc1(n1,1)=fc1(n1,1)+con2*2.
      if(n2.ge.n1)fc1(n1,n2-n1+1)=fc1(n1,n2-n1+1)+con2
      if(n3.ge.n1)fc1(n1,n3-n1+1)=fc1(n1,n3-n1+1)+con2

      if(n1.ge.n2)fc1(n2,n1-n2+1)=fc1(n2,n1-n2+1)+con2
      fc1(n2,1)=fc1(n2,1)+con2*2.
      if(n3.ge.n2)fc1(n2,n3-n2+1)=fc1(n2,n3-n2+1)+con2

      if(n1.ge.n3)fc1(n3,n1-n3+1)=fc1(n3,n1-n3+1)+con2
      if(n2.ge.n3)fc1(n3,n2-n3+1)=fc1(n3,n2-n3+1)+con2
      fc1(n3,1)=fc1(n3,1)+con2*2.

  100 continue

*Introduction of heat transfer

*        if(ntrans.ne.0)then
*        do 1200 i=1,ntrans
*        ip1=jtt(i,1)
*        ip2=jtt(i,2)
*        fk1(ip1,1)=fk1(ip1,1)+transcoef*flen3(i)/3.
*        fk1(ip2,1)=fk1(ip2,1)+transcoef*flen3(i)/3.
*          if(ip2.ge.ip1)then
*          fk1(ip1,ip2-ip1+1)=fk1(ip1,ip2-ip1+1)+transcoef*flen3(i)/6.
*          else
*          fk1(ip2,ip1-ip2+1)=fk1(ip2,ip1-ip2+1)+transcoef*flen3(i)/6.
*          endif
* 1200   continue
*        endif



*      do 6700 i=1,6
*      write(6,9990)i,(fk1(i,j),j=1,6)
* 6700 continue

*      do 6800 i=1,6
*      write(6,9990)i,(fc1(i,j),j=1,6)
* 6800 continue
* 9990 format(i5,6g12.4)

      return
      end

***********************************************************************

      subroutine tanalys(fk1,ff1,npoin,jpt,pt,npresst,fnodetemp,jbw2)
      COMMON/C1/KKK(8),NOD(9999,3),X(9998,2),LOAD(19996),
     >          JP(999),HENI(999)
      real load
      dimension fk1(9998,199),ff1(9998),jpt(999),pt(999,9980),
     >          fnodetemp(9998)

      do 400 i=1,npoin
      ff1(i)=0.
  400 continue

      do 500 i=1,npresst
      jjpt=jpt(i)
      ppt=pt(i,1)
      do 600 j=1,npoin
      jjjpt=jjpt-j+1
*        if(jjjpt.le.0)goto600
        if(jjjpt.ge.1)then
        jj=j
        jjj=jjjpt
        else
        jj=jjpt
        jjj=j-jjpt+1
        endif
*      write(*,*)'jj,jjj',jj,jjj
        if(jjj.le.jbw2)then
        ffk1=fk1(jj,jjj)
        ff1(j)=ff1(j)-ffk1*ppt
        endif
  600 continue
      fk1(jjpt,1)=1.e20
  500 continue

*      do 700 i=1,6
*      write(6,9990)ff1(i),(fk1(i,j),j=1,6)
*  700 continue
* 9990 format(1h ,7g12.4)

      call solver2(fk1,ff1,fnodetemp,npoin,jbw2,1.)
*      CALL SWEEP(fk1,ff1,fnodetemp,9998,Npoin)

      do 100 i=1,npresst
      jjpt=jpt(i)
      fnodetemp(jjpt)=fnodetemp(jjpt)+pt(i,1)
  100 continue

*      write(6,*)' '
*      write(6,*)'node, x, y, temperature'
      do 8999 i=1,npoin
*      write(6,9001)i,x(i,1),x(i,2),fnodetemp(i)
 9001 format(1h ,i5,3f10.2)
 8999 continue

*        if(istr.eq.0)then
*        do 8998 i=1,npoin
*        write(8,3237)i,x(i,1),x(i,2),fnodetemp(i),fnodepres(i)/1.e6
* 8998   continue
* 3237   FORMAT(I5,4g12.4)
*        endif

      return
      end

      SUBROUTINE SWEEP(A,B,X,N2,N)
      DIMENSION A(N2,N2),B(N2),X(N2)
      NM1=N-1
      DO 1 K=1,NM1
      KP1=K+1
      DO 2 J=KP1,N
    2 A(K,J)=A(K,J)/A(K,K)
      B(K)=B(K)/A(K,K)
      DO 3 I=KP1,N
      DO 4 J=KP1,N
    4 A(I,J)=A(I,J)-A(I,K)*A(K,J)
    3 B(I)=B(I)-A(I,K)*B(K)
    1 CONTINUE
      X(N)=B(N)/A(N,N)
      DO 5 K=1,NM1
      L=N-K
      LP1=L+1
      X(L)=B(L)
      DO 6 J=LP1,N
    6 X(L)=X(L)-A(L,J)*X(J)
    5 CONTINUE
      RETURN
      END
*************************************************
      SUBROUTINE INPUTD(ifrac,nelem,dump,tcond,npoin,jpt,pt,
     >                  npresst,expancoeff,index,reftemp,ntemp,dt,
     >                  dens,capa,cx,cy,fe2,
     >                  ntrans,transcoef,emissivity,jtt,fluidtemp,
     >                  flen3,poro,fusionheat,sun,istr,itempe,iflow,
     >                  pe1,ft1,fb1,pe2,ft2,fb2,
     >                  jpwp,pwp,ss,npwp,fpini,
     >                  accy,jbw2,ufw2,sr2,neleex,ipxc1,ipxc2,
     >                  cthre,cscale,alpha,aaalp,bbalp,ncon)

      COMMON/C1/KKK(8),NOD(9999,3),X(9998,2),LOAD(19996),
     >          JP(999),HENI(999)    
      COMMON/C2/E(9999),V(9999),THICK(9999),EARTH(9999),VOLU(9999),
     >          FOME,
     1          PHI(9999),qu(9999),RPHI(9999),Rqu(9999),TNS(9999),
     2          RTNS(9999),DIL(9999),STRF(3)
      COMMON/C3/ST(19996,998),XMK(6,6),XX(19996),BMAT(3,6,9999),
     1          STRESS(3,6,9999),XXT(19996)
      COMMON/C4/ISTAT(9999),SIGM(3,9999),STRM(3,9999),IIII,LK,FOM,KM,
     1          ORX(9999),ORY(9999),NA,NAL
	  COMMON/C6/a0(9999),c0(9999),es0(9999),akk(9999),bkk(9999),
     >          dkk(9999),tkk(9999),az0(9999),fk01(9999),alf
      DIMENSION HED(18),IZAI(9999),JPP(999),YY(999),JF(999),
     1          GE(5),GV(5),GT(5),GP(5),GC(5),GRP(5),GRC(5),      
     2          GS(5),GRS(5),GD(5),tc(5),expancoef(5),den(5),
     >          aalp(5),balp(5),aaalp(9999),bbalp(9999)
      dimension ifrac(9999),tcond(9999),jpt(999),pt(999,9980),
     >          expancoeff(9999),dens(9999),capa(9999),capac(5),
     >          cx(9999),cy(9999),fe2(9999),
     >          gfe2(5),
     >          jtt(99,2),fluidtemp(99,9980),flen3(99),poro(9999),
     >          pporo(5),sun(99,9980),
     >          fpe1(5),fft1(5),ffb1(5),
     >          fpe2(5),fft2(5),ffb2(5),
     >          pe1(9999),ft1(9999),fb1(9999),
     >          pe2(9999),ft2(9999),fb2(9999),
     >          jpwp(149),pwp(149,9980),sss(5),ss(9999),
     >          ufw1(5),ufw2(9999),sr1(5),sr2(9999)
     >          ,tcorr(999),scorr(999),
     >          alpha(9999),ga0(5),gc0(5),ges0(5),
     >          gakk(5),gbkk(5),gdkk(5),gtkk(5),gfk01(5),gec(5)


*      DOUBLE PRECISION ST
      REAL LOAD
*      character*72 hed                                                         

C   *    *    *    *    *    *    *    *    *    *    *    *    *    *  

* Jump to excavation
      IF(NA) 1,1,2

    1 WRITE(6,3000)                         
 3000 FORMAT(1h ,'STRESS  ANALYSIS  WITH  FINITE  ELEMENT METHOD')

C READING AND PRINTING OF DATA                                      

      write(6,*)
      READ (5,1001) HED
      WRITE(6,2001) HED
 1001 FORMAT(18A4)                                                      
 2001 FORMAT(1h ,20A4)

      write(6,*)
      read(5,*)istr,itempe,iflow,ielec
      write(*,*)istr,itempe,iflow,ielec

*        if((itempe.eq.2.or.itempe.eq.3).and.iflow.eq.1)then
*        write(*,*)'”ñ’èí“`”M‰ðÍ‚Æ’èíZ“§—¬‰ðÍ‚Ì˜A¬‚Í‚Å‚«‚Ü‚Ö‚ñ'
*        stop
*        endif

        if(itempe.eq.1.and.iflow.eq.2)then
        write(*,*)'’èí“`”M‰ðÍ‚Æ”ñ’èíZ“§—¬‰ðÍ‚Ì˜A¬‚Í‚Å‚«‚Ü‚Ö‚ñ'
        stop
        endif

        if((itempe.eq.2.and.iflow.eq.3).or.
     >     (itempe.eq.3.and.iflow.eq.2))then
        write(*,*)'‹«ŠEðŒˆê’è‚Ì”ñ’èí‚Æ‹«ŠEðŒ•Ï“®‚Ì”ñ’èí‚Ì˜A¬‚Í'
        write(*,*)'ƒf[ƒ^“ü—Í‚Ì“s‡ãÝ’è‚³‚ê‚Ä‚Ü‚Ö‚ñ'
        stop
        endif

        if(istr.eq.0.and.itempe.eq.0.and.iflow.eq.0)then
        write(*,*)'‚È‚É‚µ‚½‚¢‚ñ‚Å‚Á‚©H'
        stop
        endif

        if(istr.eq.0)write(6,*)'No stress analysis'
        if(istr.eq.1)write(6,*)'Elastic stress analysis'
        if(istr.eq.2)write(6,*)'Elasto-plastic stress analysis'

        if(itempe.eq.0)write(6,*)'No heat solver'
        if(itempe.eq.1)write(6,*)'Steady state heat analysis'
        if(itempe.eq.2)write(6,*)'Transient heat analysis'
        if(itempe.eq.3)write(6,*)
     >  'Transient heat analysis (variable boundary conditions)'
        if(iflow.eq.0)write(6,*)'No fluid flow'
        if(iflow.eq.1)write(6,*)'Steady fluid flow'
        if(iflow.eq.2)write(6,*)'Transient fluid flow'
        if(iflow.eq.3)write(6,*)
     >  'Transient fluid flow (variable boundary conditions)'

      write(6,*)
      READ(5,*) NPOIN,NELEM,NPRESX,NPRESY,NFORCX,NFORCY,INDEX,LYOUNG,
     >          npresst,ntrans,npwp,ipxc1,ipxc2,scalenodx,scalenody
      WRITE (6,*)'   NPOINTS     NELEM   PRESC-X   PRESC-Y'
      WRITE(6,2002)NPOIN,NELEM,NPRESX,NPRESY
      WRITE (6,*)'   NFORCEX   NFORCEY   npresst    ntrans'
      WRITE(6,2002)NFORCX,NFORCY,npresst,ntrans
      WRITE (6,*)'      Npwp'
      WRITE(6,2002)npwp
      write(6,*)'ipxc1, ipxc2',ipxc1,ipxc2

      cthre=0.
      cscale=10.

*        if(ipxc1.ne.0)then
*        write(*,*)'Input cthre, cscale. Ex. 0., 1.'
*        write(*,*)'Increase cscale to, ex. 10., if vibration is observed.'
*        read(*,*)cthre,cscale
*        endif

* 1002 FORMAT(10I5)
 2002 FORMAT(1h ,4I10)

      write(6,*)
      IF(INDEX.EQ.0) GO TO 2500                                         
      WRITE(6,2401)                                                     
 2401 FORMAT(1h ,'PLANE STRESS CONDITION')
      GO TO 2550                                                        
 2500 WRITE(6,2402)                                                     
 2402 FORMAT(1h ,'PLANE STRAIN CONDITION')

 2550 write(6,*)
      WRITE (6,2403) LYOUNG                                             
 2403 FORMAT(1h ,'  NUMBER OF MATERIALS ='I3)                           
      IW=2*NPOIN

      write(6,*)
      WRITE(6,3003)                                                     
 3003 FORMAT(1h ,'  NELEM NOD(I,1) NOD(I,2) NOD(I,3) NMAT')   
      READ(5,*)((NOD(I,J),J=1,3),IZAI(I),I=1,NELEM)

 1015 FORMAT(16I4)
      do 7266 i=1,nelem
      WRITE(6,2004)I,(NOD(I,J),J=1,3),IZAI(I)
 2004 FORMAT(1h ,I7,2X,3I4,2X,I2)                                        
 7266 continue
 1004 FORMAT(15I5)                                                      

      write(6,*)
      WRITE(6,3102)
 3102 FORMAT(1h ,'NPOIN X-CO-OD(I) Y-CO-OD(I)')
      READ (5,*) ((X(I,J),J=1,2),I=1,NPOIN)

      do 2688 i=1,npoin
*      do 2688 j=1,2
      x(i,1)=x(i,1)*scalenodx
      x(i,2)=x(i,2)*scalenody
 2688 continue

      do 7267 i=1,npoin
      WRITE(6,2003)I,(X(I,J),J=1,2)
 7267 continue
 1003 FORMAT(8F8.3)                                                     
 2003 FORMAT(1h ,I5,1X,2F9.4,2X)

      do 8701 i=1,nelem
      sumx=0.
      sumy=0.
      do 8702 j=1,3
      sumx=sumx+x(nod(i,j),1)
      sumy=sumy+x(nod(i,j),2)
 8702 continue
      cx(i)=sumx/3.
      cy(i)=sumy/3.
 8701 continue

      write(6,*)' '
      write(6,*)'Elements'' center'
      do 8704 i=1,nelem
      write(6,8703)i,cx(i),cy(i)
 8703 format(1h ,i5,2g12.4)
 8704 continue

C BOUNDARY CONDITION                                                
C PRESCRIBED DISPLACEMENT                                           
      NPRESC=NPRESX+NPRESY                                              
*      WRITE(6,2410)                                                     
* 2410 FORMAT(1h ,' *BOUNDARY CONDITIONS*')                              
*      WRITE(6,2411)                                                     
*      WRITE(6,2412)                                                     
 2411 FORMAT(1h ,'  *PRESCRIBED DISPLACEMENT')                          
 2412 FORMAT(1H ,19H   - X-DIRECTION - )                                
 2413 FORMAT(1H ,16H    (NOD, DISP ))                                   

      write(6,*)
      WRITE(6,*)'Prescribed displacement node, disp'
        IF(NPRESX.EQ.0) GO TO 73
*      write(*,*)npresx
      READ(5,*) (JPP(I),I=1,NPRESX)
      READ(5,*) (YY(I),I=1,NPRESX)
      DO 50 I=1,NPRESX
      JP(I)=JPP(I)*2-1                                                  
   50 HENI(I)=YY(I)                                                     

      write(6,*)'X-direction'
      do 7241 i=1,npresx
      WRITE(6,2415) JPP(I),HENI(I)
 2415 FORMAT(I8,F8.3)
 7241 continue

   73 continue


        IF(NPRESY.EQ.0) GO TO 72                                          
      WRITE(6,*)'Y-direction'
      READ(5,*) (JPP(I),I=1,NPRESY)                                  
      READ(5,*) (YY(I),I=1,NPRESY)                                   
      DO 51 I=1,NPRESY                                                  
      K=I+NPRESX                                                        
      JP(K)=JPP(I)*2                                                    
   51 HENI(K)=YY(I)
      do 7543 i=1,npresy
      WRITE(6,2415) JPP(I),HENI(NPRESX+I)
 7543 continue

C *** FORCE-TYPE BOUNDARY CONDITION                                     
   72 DO 45  I=1,IW                                                     
   45 LOAD(I)=0.                                                        
      NFORCE=NFORCX+NFORCY                                              

        IF(NFORCE) 400,400,802                                            
  802 continue
      write(6,*)
      WRITE(6,*)'Nodal point load'

        IF(NFORCX.EQ.0) GO TO 52                                          
      READ(5,*) (JPP(I),I=1,NFORCX)                                  
      READ(5,*) (YY(I),I=1,NFORCX)                                   
      DO 53 I=1,NFORCX                                                  
      JF(I)=2*JPP(I)-1                                                  
   53 LOAD(JF(I))=YY(I)                                                 

      WRITE(6,*)'X-direction'
      do 7544 i=1,nforcx
      WRITE(6,2422) JPP(I),YY(I)
 7544 continue

 2422 FORMAT(1H ,i5,g12.4)

   52   IF(NFORCY.EQ.0) GO TO 59                                          
      READ(5,*) (JPP(I),I=1,NFORCY)                                  
      READ(5,*) (YY(I),I=1,NFORCY)                                   
 1013 FORMAT(4F16.3)                                                    
      DO 54 I=1,NFORCY                                                  
      K=NFORCX+I                                                        
      JF(K)=2*JPP(I)                                                    
   54 LOAD(JF(K))=YY(I)                                                 

      WRITE(6,*)'Y-direction'
      do 7545 i=1,nforcy
      WRITE(6,2422) JPP(I),YY(I)
 7545 continue

   59 CONTINUE                                                          
  400 CONTINUE                                                          


C                                                                       
C *** MATERIAL PROPERTIES                                               
C 

 
      DO 70 I=1,LYOUNG
      read(5,*) ga0(I),gc0(I),ges0(I),gakk(I),gbkk(I),gdkk(I),
     >          gtkk(I),gfk01(I),gec(I),ncon
 1080 FORMAT(1H ,8g10.3)
      READ(5,*) GE(I),GV(I),GT(I),GP(I),GC(I),GRP(I),GRC(I),
     1          GS(I),GRS(I),GD(I),
     >          tc(i),capac(i),expancoef(i),den(i),
     >          gfe2(i),pporo(i),
     >          fpe1(i),fft1(i),ffb1(i),fpe2(i),fft2(i),ffb2(i),
     >          sss(i),
     >          ufw1(i),sr1(i),aalp(i),balp(i)

      write(6,*) GE(I),GV(I),GT(I),GP(I),GC(I),GRP(I),GRC(I),
     1          GS(I),GRS(I),GD(I),
     >          tc(i),capac(i),expancoef(i),den(i),
     >          gfe2(i),pporo(i),
     >          fpe1(i),fft1(i),ffb1(i),fpe2(i),fft2(i),ffb2(i),
     >          sss(i),
     >          ufw1(i),sr1(i),aalp(i),balp(i)

      write(6,*)
      WRITE(6,*)'Material',I
      write(6,*)
      WRITE(6,7293)
 7293 FORMAT(1h ,'         E         V     THICK',
     >           '       PHI       UCS      RPHI')
      write(6,7291)GE(I),GV(I),GT(I),GP(I),GC(I),GRP(I)
 7291 FORMAT(1H ,6g10.3)
      write(6,*)
      write(6,7294)
 7294 format(1h ,'     R-UCS       T0       rt0       dil',
     >           '     th-coSpec-heat  expancof       den')
      write(6,7292)GRC(I),GS(I),GRS(I),GD(I),tc(i),capac(i),
     >             expancoef(i),den(i)
 7292 format(1h ,8g10.3)

      write(6,*)
      write(6,7295)gfe2(i),pporo(i),fpe1(i),fft1(i),ffb1(i),fpe2(i),
     >             fft2(i),ffb2(i),
     >             sss(i),ufw1(i),sr1(i),aalp(i),balp(i)
 7295 format(1h ,'Freeze completes at       ',f10.2,/,
     >       1h ,'Effective porosity        ',f10.5,/,
     >       1h ,'Pearmeability at T1 (m^2) ',e12.4,/,
     >       1h ,'T1 (oC)                   ',f10.2,/,
     >       1h ,'b1 (/MPa)                 ',e12.4,/,
     >       1h ,'Pearmeability at T2 (m^2) ',e12.4,/,
     >       1h ,'T2 (oC)                   ',f10.2,/,
     >       1h ,'b2 (/MPa)                 ',e12.4,/,
     >       1h ,'Specific storage (m^-1)   ',e12.4,/,
     >       1h ,'Unfrozen water ratio      ',f10.4,/,
     >       1h ,'Saturation ratio          ',f10.4,/,
     >       1h ,'aalp                      ',f10.4,/,
     >       1h ,'balp (MPa^-1)             ',f10.4,/)
   70 CONTINUE
C                                                                       
C     INITIAL STRESSES AND LIMIT OF ITERATION                           
C                                                                       

      READ(5,*) NLIM,dump,fome,NAL,npoiex,neleex,npresxex
      write(6,8171)NLIM,dump,fome,NAL,npoiex,neleex,npresxex
 8171 format(1h ,'NLIM,dump,fome,NAL,npoiex,neleex,npresxex',
     >       7g10.2)
      READ(5,*) (STRF(I),I=1,3),accy
      write(6,8172)(strf(i),i=1,3),accy
 8172 format(1h ,'strf, accy',4g12.4)

      DO 71 K=1,NELEM
      I=IZAI(K)
      ISTAT(K)=I
	  a0(k)=ga0(I)
	  c0(k)=gc0(I)
	  es0(k)=ges0(I)
	  akk(k)=gakk(I)
	  bkk(k)=gbkk(I)
	  dkk(k)=gdkk(I)
	  tkk(k)=gtkk(I)
	  az0(k)=gec(I)
	  fk01(k)=gfk01(I)
      E(K)=GE(I)
      V(K)=GV(I)
      THICK(K)=GT(I)
      EARTH(K)=accy*den(i)
      PHI(K)=GP(I)                                                      
      qu(K)=GC(I)                                                      
      RPHI(K)=GRP(I)                                                    
      Rqu(K)=GRC(I)                                                    
      TNS(K)=GS(I)                                                      
      RTNS(K)=GRS(I)                                                    
      DIL(K)=GD(I)
      ifrac(k)=0
      capa(k)=capac(i)
      tcond(k)=tc(i)
      dens(k)=den(i)
      expancoeff(k)=expancoef(i)
      fe2(k)=gfe2(i)
      poro(k)=pporo(i)

      pe1(k)=fpe1(i)
      ft1(k)=fft1(i)
      fb1(k)=ffb1(i)
      pe2(k)=fpe2(i)
      ft2(k)=fft2(i)
      fb2(k)=ffb2(i)

      ss(k)=sss(i)
      ufw2(k)=ufw1(i)
      sr2(k)=sr1(i)
      aaalp(k)=aalp(i)
      bbalp(k)=balp(i)
* Initial value
      alpha(k)=aaalp(k)
   71 CONTINUE                                                          

      write(6,*)
      WRITE(6,3301) (STRF(I),I=1,3)
 3301 FORMAT(1h ,'INITIAL STRESSES',3g12.3)
      write(6,*)
      write(6,6221)accy
 6221 format(1h ,'Acceleration in -y directon        ',f12.4)

      write(6,*)
      WRITE(6,3302) NAL,npoiex,neleex,NLIM,FOME
 3302 FORMAT(1h ,'NUMBER OF MINING STEPS =                ',I2,/,
     >       1h ,'Number of points excavated              ',I5,/,
     >       1h ,'Number of elements excavated            ',I5,/,
     >       1h ,'MAXIMUM NUMBER OF ITERATIONS =          ',I6,/,
     >       1h ,'MAXIMUM DIFERENCE OF NODAL POINT LOADS =',g12.1)
      DO 75 I=1,NELEM
      PHI(I)=3.141593*PHI(I)/180.0
   75 RPHI(I)=3.141593*RPHI(I)/180.0
      GO TO 3
C   *    *    *    *    *    *    *    *    *    *    *    *    *    *  
C     DETERMINE BAND WIDTH                                              
*     *     *      *        *           *          *        *         * 
***********************************
*Excavation

    2 NPOIN=NPOIN-npoiex
      NELEM=NELEM-neleex
      IW=2*NPOIN
      NPRESx=NPRESx-npresxex
      NPRESC=NPRESX+NPRESY
***********************************

    3 J=0                                                               
      DO 340 N=1,NELEM                                                  
      DO 340 I=1,3                                                      
      DO 325 L=1,3                                                      
      KK=IABS(NOD(N,I)-NOD(N,L))                                        
      IF(KK-J) 325,325,320                                              
  320 J=KK                                                              
  325 CONTINUE                                                          
  340 CONTINUE                                                          
      JBW=2*J+2
      jbw2=jbw/2

      KKK(1)=NPOIN
      KKK(2)=NELEM
      KKK(3)=NPRESC
      KKK(4)=NFORCE
      KKK(5)=IW
      KKK(6)=JBW
      KKK(7)=NLIM
      KKK(8)=INDEX

        if(itempe.eq.0.and.iflow .eq.0)then
        ntemp=1
        else
        read(5,*)reftemp,dt,ntemp,fusionheat,fpini
        write(6,*)
        write(6,6901)reftemp,dt,ntemp,fusionheat,fpini
 6901   format(1h ,'Reference temperature              ',f12.4,/,
     >         1h ,'Time step for transient heat solver',g12.4,/,
     >         1h ,'Number of time steps               ',i12,/,
     >         1h ,'Fusion heat in W/m^3               ',g12.4,/,
     >         1h ,'Initial pore presure in Pa         ',g12.4)
        endif

        if(itempe.eq.0.and.iflow .eq.1)ntemp=1

* Prescribed temperature for nodes

        if(npresst.ne.0)then
          if(itempe.eq.3)then
*          do 8122 i=1,npresst
*          read(5,*)jpt(i),(pt(i,j),j=1,ntemp)
* 8122     continue
*-----------------------------------------
          read(5,*)(jpt(i),i=1,npresst)
          do 8123 j=1,ntemp
          read(5,*)pt(1,j)
*            if(j.gt.690)then
*            write(*,*)'itemp, temp',j, pt(1,j)
*            endif
          do 8123 i=2,npresst
          pt(i,j)=pt(1,j)
 8123     continue
*-----------------------------------------
          endif
          if(itempe.eq.1.or.itempe.eq.2)then
          do 8622 i=1,npresst
          read(5,*)jpt(i),pt(i,1)
          do 8622 j=1,ntemp
          pt(i,j)=pt(i,1)
 8622     continue
          endif

        write(6,*)
        write(6,*)' Node and prescribed temperature'
        do 7193 i=1,npresst
        write(6,7194)jpt(i)
        do 9239 j=1,ntemp
        write(6,9240)pt(i,j)
 9239   continue
 7193   continue
 7194   format(1h ,i5)
 9240   format(1h ,f10.2)
        endif

* Prescribed temperature for surrounding fluid

        if(ntrans.ne.0)then
        read(5,*)transcoef,emissivity
        write(6,1992)transcoef,emissivity
 1992   format(1h ,'transcoeff                          ',g12.5,/,
     >         1h ,'emissivity                          ',g12.5)

        write(5,*)' ntrans=',ntrans

          do 4622 i=1,ntrans
          read(5,*)jtt(i,1),jtt(i,2),tcorr(i),scorr(i)
 4622     continue

*          write(*,*)'itempe,ntemp',itempe,ntemp

          do 4621 j=1,ntemp
            if(itempe.eq.2.and.j.ge.2)goto4620
          read(5,*)fluidtemptemp,suntemp
          write(*,*)'fluidtemptemp,suntemp',fluidtemptemp,suntemp
          do 4621 i=1,ntrans
          fluidtemp(i,j)=fluidtemptemp*tcorr(i)
          sun(i,j)=suntemp*scorr(i)
 4621     continue

 4620     continue
*              write(*,*)'OOk'

        write(6,*)
        write(6,*)' Node and Prescribed temperature',
     >            ' for sorrounding fluid'
        do 2193 i=1,ntrans
        write(6,2194)jtt(i,1),jtt(i,2)
        xx1=x(jtt(i,1),1)
        yy1=x(jtt(i,1),2)
        xx2=x(jtt(i,2),1)
        yy2=x(jtt(i,2),2)
        flen3(i)=sqrt((xx2-xx1)**2+(yy2-yy1)**2)
        do 2239 j=1,ntemp
        write(6,9240)fluidtemp(i,j)
 2239   continue
 2193   continue
 2194   format(1h ,2i5)
        endif

* Density and viscosity are automatically calculated
* for the temperature
* assuming the fluid is water and temperature.
*        if(iflow.ne.0)then
*        read(5,*)fdens,fmyu
*        write(6,4848)fdens,fmyu
* 4848   format(1h ,'Density of pore fluid              ',g12.4,/,
*     >         1h ,'Viscosity of pore fluid            ',g12.4)
*        endif

* Prescribed pore fluid pressure for nodes

        if(npwp.ne.0)then

          if(iflow.eq.3)then
          do 1122 i=1,npwp
          read(5,*)jpwp(i),(pwp(i,j),j=1,ntemp)
 1122     continue
          endif

          if(iflow.eq.1.or.iflow.eq.2)then
          do 1622 i=1,npwp
          read(5,*)jpwp(i),pwp(i,1)
          do 1622 j=1,ntemp
          pwp(i,j)=pwp(i,1)
 1622     continue
          endif

        write(6,*)
        write(6,*)' Node and prescribed pore fluid pressure'
        do 1193 i=1,npwp
        write(6,1194)jpwp(i)
        do 1239 j=1,ntemp
        write(6,1240)pwp(i,j)
 1239   continue
 1193   continue
 1194   format(1h ,i5)
 1240   format(1h ,e12.4)

        endif

      write(6,*)
      WRITE(6,3021) JBW                                                 
 3021 FORMAT(1h ,'BAND WIDTH =',I6,' IS ALLOWABLE WITHIN 998')        
      IF(998-JBW) 299,311,311                                            
  299 WRITE(6,3011) JBW                                                 
 3011 FORMAT(1h ,'BAND WIDTH=',I6,' EXCEEDS ALLOWABLE NUMBER 998')    
  300 STOP                                                              
  311 CONTINUE                                                          

      RETURN                                                            
      END                                                               
*****************************************************************
      SUBROUTINE STIFNS(c1,itemp,icon,prstr1,prstr2,prstr3)
      COMMON/C1/KKK(8),NOD(9999,3),X(9998,2),LOAD(19996),
     >          JP(999),HENI(999)    
      COMMON/C2/E(9999),V(9999),THICK(9999),EARTH(9999),VOLU(9999),FOME,     
     1          PHI(9999),qu(9999),RPHI(9999),Rqu(9999),TNS(9999),         
     2          RTNS(9999),DIL(9999),STRF(3)                              
      COMMON/C3/ST(19996,998),XMK(6,6),XX(19996),BMAT(3,6,9999),              
     1          STRESS(3,6,9999),XXT(19996)                                
      COMMON/C4/ISTAT(9999),SIGM(3,9999),STRM(3,9999),IIII,LK,FOM,KM,      
     1          ORX(9999),ORY(9999),NA,NAL                                
      DIMENSION XE(3,2),c1(3,3,9999),
     >          prstr1(9999),prstr2(9999),prstr3(9999)
*      DOUBLE PRECISION ST
      REAL LOAD                                                         

      IW=KKK(5)                                                         
      JBW=KKK(6)                                                        
      NELEM=KKK(2)                                                      
      DO 60 I=1,IW                                                      
      DO 60 J=1,JBW
   60 ST(I,J)=0.                                                        
C   *    *    *    *    *    *    *    *    *    *    *    *    *    *  
   93 DO 80 LK=1,NELEM                                                  
      DO 85 I=1,3                                                       
      JJ=NOD(LK,I)                                                      
      XE(I,1)=X(JJ,1)                                                   
   85 XE(I,2)=X(JJ,2)                                                   
C   *    *    *    *    *    *    *    *    *    *    *    *    *    *  
C     CALCULATION OF ELEMENT STIFFNESS AND STRESS MATRICES              
C   *    *    *    *    *    *    *    *    *    *    *    *    *    *  
      CALL FEM (XE,c1,itemp,icon,prstr1,prstr2,prstr3)                                                     
      DO 5 LL=1,3                                                       
      DO 5 KK=1,3                                                       
      M=2*(NOD(LK,KK)-1)                                                
      N=2*(NOD(LK,LL)-1)                                                
      I=2*(KK-1)                                                        
      J=2*(LL-1)                                                        
      DO 5 NJ=1,2                                                       
      DO 5 MI=1,2                                                       
      MMI=M+MI                                                          
      NNJ=N+NJ-MMI+1                                                    
      IF(NNJ.LE.0) GO TO 5                                              
      IMI=I+MI                                                          
      JNJ=J+NJ                                                          
      ST(MMI,NNJ)=ST(MMI,NNJ)+XMK(IMI,JNJ)                              
    5 CONTINUE                                                          
   80 CONTINUE                                                          
      RETURN                                                            
      END                                                               
**********************************************************
      SUBROUTINE ANALYS(ifrac,dump,elethermstrain,c1,npoin,istr,
     >                  neleex,
     >                  ipxc1,ipxc2,
     >                  cthre,cscale,stsv,elepres,alpha,
     >                  pload,tload,prstr1,prstr2,prstr3)
      COMMON/C1/KKK(8),NOD(9999,3),X(9998,2),LOAD(19996),
     >          JP(999),HENI(999)    
      COMMON/C2/E(9999),V(9999),THICK(9999),EARTH(9999),VOLU(9999),FOME,     
     1          PHI(9999),qu(9999),RPHI(9999),Rqu(9999),TNS(9999),         
     2          RTNS(9999),DIL(9999),STRF(3)                              
      COMMON/C3/ST(19996,998),XMK(6,6),XX(19996),BMAT(3,6,9999),              
     1          STRESS(3,6,9999),XXT(19996)                                
      COMMON/C4/ISTAT(9999),SIGM(3,9999),STRM(3,9999),IIII,LK,FOM,KM,      
     1          ORX(9999),ORY(9999),NA,NAL                                
      COMMON/C5/QE(2,9999),QS(6,9999),PRST(3,9999),NC
*	  COMMON/C6/a0,c0,es0,d0,s0,al0,az0,fk01,z
      dimension ifrac(9999),elethermstrain(9999),c1(3,3,9999),
     >          strmr(3,9999),
     >          prstr1(9999),prstr2(9999),prstr3(9999),
     >          str(19996,998),stsv(9999),elepres(9999),
     >          alpha(9999),pload(19996),tload(19996)

*      DOUBLE PRECISION ST
      REAL LOAD                                                         
C                                                                       

      NPOIN=KKK(1)                                                      
      NELEM=KKK(2)                                                      
      NLIM=KKK(7)                                                       
      IW=KKK(5)                                                         
      IF(NA) 9,9,99                                                     
C     INITIAL CONDITION                                                 
    9 DO 1 I=1,NELEM                                                    
      DO 2 J=1,3                                                        
      SIGM(J,I)=STRF(J)
      STRM(J,I)=0.0
    2 CONTINUE

    1 CONTINUE                                                          
      DO 4 I=1,IW                                                       
    4 XXT(I)=0.0                                                        
   99 IIII=0                                                            
C     NC=0

*Iteration for elasto-plastic

 1000 CONTINUE                                                          

    7 continue
      iw=kkk(5)
      jbw=kkk(6)
      do 6444 i=1,iw

*      xxt(i)=0.

      do 6444 j=1,jbw
      str(i,j)=st(i,j)
 6444 continue

      do 6443 i=1,iw
*      load(i)=rload(i)+pload(i)+tload(i)

*      write(*,*)i,load(i)

      load(i)=load(i)+pload(i)+tload(i)
 6443 continue

      CALL SOLVER

*121--------------------------------------------------------
      if(ipxc1.ne.0)then
        du1=xxt(2*ipxc1-1)-xxt(2*ipxc2-1)
        write(*,*)' du1',du1
        if(du1.gt.cthre)then
*        if(du1.gt.1.e-6)then
*        if(du1.gt.1.e-5)then
*        if(du1.gt.1.e-4)then
        write(*,*)' du1.gt.0!'
*        DO 6450 I=1,NELEM
*        DO 6450 J=1,3
*        SIGM(J,I)=STRF(J)
*        STRM(J,I)=0.0
* 6450   CONTINUE
        do 6447 i=1,iw
*        xxt(i)=0.
        load(i)=0.
 6447   continue
*        cscale=1.
        load(ipxc1*2-1)=cscale
        load(ipxc2*2-1)=-cscale
        do 6445 i=1,iw
        do 6445 j=1,jbw
        st(i,j)=str(i,j)
 6445   continue
*        iiii=1
        CALL solver
        du2=xxt(2*ipxc1-1)-xxt(2*ipxc2-1)
        fpxc=-du1/(du2-du1)*cscale
*        fpxc=-du1/du2*cscale
        write(*,*)'du1,du2,fpxc',du1,du2,fpxc

*        DO 6449 I=1,NELEM                                                    
*        DO 6449 J=1,3                                                        
*        SIGM(J,I)=STRF(J)
*        STRM(J,I)=0.0
* 6449   CONTINUE
        do 6448 i=1,iw
*        xxt(i)=0.
        load(i)=0.
 6448   continue
        load(ipxc1*2-1)=fpxc
        load(ipxc2*2-1)=-fpxc
        do 6446 i=1,iw
        do 6446 j=1,jbw
        st(i,j)=str(i,j)
 6446   continue
*        iiii=1
        call solver
        endif
      endif
*121------------------------------------------------------------
*        if(ipxc1.ne.0.and.ipxc.ne.2)return

*20140218
*      DO 3 I=1,IW
*      LOAD(I)=0.0
*    3 continue
*20140218

      FOM=0.0
      FON=0.0
    6 continue

*        if(iiii.eq.0)then
*          write(6,*)
*          WRITE(6,3389)
* 3389     FORMAT(1h ,'OUTPUT AT first STEP')
*          write(6,*)
*          WRITE(6,3398)
*          endif
   56 DO 10 LK=1,NELEM
*     do 2022 k444=1,20                                     
      CALL ANAROK(ifrac,dump,elethermstrain,c1,strmr,istr,
     >            prstr1,prstr2,prstr3,ipxc1,stsv,elepres,
     >            alpha,itemp)
* 2022 continue
      IF( FOM.GT.FON )FON=FOM
   10 CONTINUE
C
C     KAI-NO SHU-SOKU HANTEI

      write(*,6161)iiii+1,fon
 6161 format(1h ,'iiii+1, fon',i6,g12.4)

      if(iiii+1.lt.nlim)then
        IF( FON.GT.FOME ) GO TO 50
      endif

      write(6,*)' '
      WRITE(6,201) IIII+1
  201 FORMAT(1h ,'SOLUTION AT ITERATION =',I4)
   11 continue

*      write(*,*)'elemout in analys'
*      call elemout(cx,cy,strmr,prstr1,prstr2,prstr3,
*     >             ifrac,flowv,nelem)
      write(6,*)' '
      WRITE(6,3487)

      IF(NA-NAL) 70,71,71
   70 NA=NA+1
      DO 20 I=1,IW
   20 LOAD(I)=0.0

*      DO 21 LK=NELEM-1,NELEM

      do 21 lk=nelem-neleex,nelem
      CALL ANACOL
   21 continue

      WRITE(6,3388) NA
 3388 FORMAT(/////'   * FACE ADVANCING 1 STEP *'/'     NUMBER OF MINING 
     1STEPS ='I3)
      GO TO 72
   71 NA=0
   72 RETURN
   50 continue
      IIII=IIII+1
      IF( IIII.LT.NLIM ) GO TO 30

      WRITE(6,202) NLIM
  202 FORMAT(//'     SOLUTION after ITERATION ='I4)

   13 continue

      WRITE(6,3487)
 3487 FORMAT(1h ,'Nodal displacement')
      do 3644 i=1,npoin
      WRITE(6,3497)I,(XXT(2*(I-1)+J),J=1,2)
 3644 continue

 3497 FORMAT(6(I5,2f12.6))
      STOP                                                              
   30 CONTINUE                                                          
      GO TO 1000
*      stop
      END
*********************************************************
      SUBROUTINE ANAROK(ifrac,dump,elethermstrain,c1,strmr,istr,
     >                  prstr1,prstr2,prstr3,ipxc1,stsv,elepres,
     >                  alpha,itemp)
      COMMON/C1/KKK(8),NOD(9999,3),X(9998,2),LOAD(19996),
     >          JP(999),HENI(999)    
      COMMON/C2/E(9999),V(9999),THICK(9999),EARTH(9999),VOLU(9999),FOME,     
     1          PHI(9999),qu(9999),RPHI(9999),Rqu(9999),TNS(9999),         
     2          RTNS(9999),DIL(9999),STRF(3)                              
      COMMON/C3/ST(19996,998),XMK(6,6),XX(19996),BMAT(3,6,9999),              
     1          STRESS(3,6,9999),XXT(19996)                                
      COMMON/C4/ISTAT(9999),SIGM(3,9999),STRM(3,9999),IIII,LK,FOM,KM,      
     1          ORX(9999),ORY(9999),NA,NAL                                
      COMMON/C5/QE(2,9999),QS(6,9999),PRST(3,9999),NC 
      COMMON/C6/a0(9999),c0(9999),es0(9999),akk(9999),bkk(9999),
     >          dkk(9999),tkk(9999),az0(9999),fk01(9999),alf  
      DIMENSION SIGT(3),IIL(6),DIS(6),c1(3,3,9999)
      dimension ifrac(9999),elethermstrain(9999),
     >          strmr(3,9999),prstr1(9999),prstr2(9999),prstr3(9999),
     >          stsv(9999),elepres(9999),alpha(9999),sr(6,9999)
*
*      DOUBLE PRECISION ST
      REAL LOAD                                                         
      pai=3.141593

*      write(*,*)'Entering anarok'
                                                                       
      DO 7 I=1,3                                                        
      IIL(2*I-1)=2*NOD(LK,I)-1                                          
    7 IIL(2*I)=2*NOD(LK,I)

      DO 182 I=1,6
      K=IIL(I)
      DIS(I)=XXt(K)
  182 continue

*      write(*,*)'lk, dis',lk,(dis(i),i=1,6)

*      write(*,*)'dis',(dis(i),i=1,6)

      DO 184 I=1,3
*’e‘Y«‚É‚Í‚±‚±‚Å‰Šú‰»‚·‚é•K—v‚ª‚ ‚éB

*      write(*,*)'iiii',iiii

*‚Ð‚¸‚Ý‚ª‚¿‚á‚ñ‚ÆŒvŽZ‚³‚ê‚Ä‚È‚¢‚ñ‚Å‚ÍH

*      DO 17 I=1,3
*        if(iiii.eq.0)then
*        sigt(i)=0.
*        else
*        SIGT(I)=SIGM(I,LK)
*        endif
*   17 continue

        if(iiii.eq.0.and.ipxc1.eq.0)strmr(i,lk)=0.
      strmr(i,lk)=strmr(i,lk)+strm(i,lk)
      strm(i,lk)=0.
      DO 184 J=1,6
      STRM(I,LK)=STRM(I,LK)+BMAT(I,J,LK)*DIS(J)
  184 continue

*      write(*,*)'lk,strm1,strm2,strm3 at 2697',
*     >           lk,strm(1,lk),strm(2,lk),strm(3,lk)

* strm is correct up to here.

      epsx=strm(1,lk)
      epsy=strm(2,lk)
      epsxy=strm(3,lk)/2.

      call sig(epsx,epsy,epsxy,eps1,eps2,alf)

*      write(*,*)'eps1,eps2',lk,eps1,eps2
      prstr1(lk)=eps1
      prstr2(lk)=eps2
      prstr3(lk)=alf

*      epsx=strm(1,lk)+strmr(1,lk)
*      epsy=strm(2,lk)+strmr(2,lk)
*      epsxy=(strm(3,lk)+strmr(3,lk))/2.
*  	  fk0(lk)=(d0/pai0*(atan((SIGM(2,lk)-s0)/aL0)
*    >   +pai0/2.)+az0)

*      write(*,*)'lk,strm1,strm2,strm3 at 2800',
*     >           lk,strm(1,lk),strm(2,lk),strm(3,lk)

*       write(*,*)'sigm before',sigm(1,lk),sigm(2,lk)

*      write(*,*)'ppstr',ppstr(lk)
*        if(itemp.ne.0)then
        do 183 i=1,3
*          write(6,*)'iiii,lk',iiii,lk
*          if(iiii.eq.0)then
          sigm(i,lk)=0.
*          endif

*      write(*,*)'lk,strm,elethermstrain',lk,(strm(j,lk),j=1,2),
*     >          elethermstrain(lk)

        DO 183 j=1,3
*          if(iiii.eq.0)then
            if(j.eq.1.or.j.eq.2)then
            SIGM(I,LK)=SIGM(I,LK)+c1(i,j,lk)
*     >                *(strm(j,lk)-elethermstrain(lk)-ppstr(lk))
     >                *(strm(j,lk)-elethermstrain(lk))
*     >                *strm(j,lk)
            else

*        write(*,*)'i,lk,j,sigm,c1',i,lk,j,sigm(i,lk),c1(i,j,lk)

            SIGM(I,LK)=SIGM(I,LK)+c1(i,j,lk)*strm(j,lk)
            endif
*          endif
  183   continue
*      write(6,*)'lk,SIGM(2,LK)',lk,SIGM(2,LK),strm(2,LK),ad02(lk)
*      write(6,*)'lk,SIGM(1,LK)',lk,SIGM(1,LK),strm(1,LK),ad01(lk)
*       E1=SIGM(2,LK)/strm(2,lk)   
*       E(lk)=E1
*       write(6,*)'lk,E(LK),V1',lk,E(LK),V1,SIGM(2,LK),strm(2,LK)
*      sigm(1,lk)=sigm(1,lk)+alpha(lk)*elepres(lk)
*      sigm(2,lk)=sigm(2,lk)+alpha(lk)*elepres(lk)

*       write(*,*)'sigm after',sigm(1,lk),sigm(2,lk)

*        else
*‚â‚¯‚­‚»‚¾*************
*      DO 9183 I=1,3
*      DO 9183 J=1,6
*      SIGM(I,LK)=SIGM(I,LK)+STRESS(I,J,LK)*DIS(J)
* 9183 continue
*************************
*        endif


*      write(6,*)(strm(i,lk),i=1,3) 

      sx =sigm(1,lk)
      sy =sigm(2,lk)
      sxy=sigm(3,lk)

      do 2645 i=1,3
      SIGT(I)=SIGM(I,LK)
 2645 continue
*--------------------------
*      write(*,9333)lk,sx,sy,sxy
* 9333 format(1h ,'lk,sx,sy,sxy',i5,3g12.4)
*---------------------------

      call sig(sx,sy,sxy,ps1,ps2,alf)

*        IF(IIII.eq.0)then
*        WRITE(6,3499) iiii,LK,(STRM(I,LK)*1.e6,I=1,3),
*     >               (SIGM(I,LK)/1.e6,I=1,3),
*     >                PRST(1,LK)/1.e6,prst(2,lk)/1.e6,prst(3,lk),
*     >                Ifrac(LK)
* 3499 FORMAT(1H ,2I4,3f10.0,5f10.2,f7.2,I5)
*        endif

*      write(6,*)lk,ps1,ps2
*
**********************************

*        if(itemp.eq.1)then
*        ps11=ps1+elepres(lk)*alpha(lk)
*        ps22=ps2+elepres(lk)*alpha(lk)
*        endif

      ps11=ps1
      ps22=ps2
*      write(*,*)'ps1,ps2,elepres',ps1,ps2,elepres(lk)

* Modified stress severity
      sst=-ps11/tns(lk)
      cc =(1.+sin(phi(lk)))/(1.-sin(phi(lk)))
      ssc=(ps22-ps11)/(-qu(lk)+ps11*cc-ps11)
        if(ssc.gt.-sst)then
        stsv(lk)=ssc
        else
        stsv(lk)=sst
        endif
**********************************
        if(istr.ne.2)then
*          if(itemp.eq.1)then
*          PRST(1,LK)=PS1-elepres(lk)*alpha(lk)
*          PRST(2,LK)=PS2-elepres(lk)*alpha(lk)
*          else
          PRST(1,LK)=PS1
          PRST(2,LK)=PS2
*          endif

        PRST(3,LK)=alf
       E1=PRST(2,LK)/prstr2(lk)
       V1=-prstr1(lk)/prstr2(lk)
       V(lk)=V1	   
       E(lk)=E1
*      write(6,*)'lk,E1',lk,E(lk),V(lk) 

*      write(*,*)'prst1,prst2',prst(1,lk),prst(2,lk)

* Strength reqired not to fail

        goto8111
        endif

* Determination of element failure and assignment of stresses
*
* ifrac   0    intact
*        -1    tensile failure
*        -2    compressive failure
*        -3    double tensile
*        -4    tensile then compressive
*

        ifracr=ifrac(lk)

* Tensile failure initiation

*      write(6,*)'lk,ps1,ps2,tns(lk),qu(lk)',
*     >           lk,ps1,ps2,tns(lk),qu(lk)

        if(ifrac(lk).eq.0.and.ps22.gt.tns(lk))then
        ifrac(lk)=-3
        endif
        if(ifrac(lk).eq.-3)then
*        ps1=rtns(lk)+alpha(lk)*elepres(lk)
*        ps2=rtns(lk)+alpha(lk)*elepres(lk)
      ps1=rtns(lk)
      ps2=rtns(lk)
        endif

        if(ifrac(lk).eq.0.and.ps11.gt.tns(lk))then
        ifrac(lk)=-1
        endif
        if(ifrac(lk).eq.-1)then
*        ps1=rtns(lk)+alpha(lk)*elepres(lk)
      ps1=rtns(lk)
        endif


* Compressive failure
*      cc =(1.+sin( phi(lk)*pai/180.))/(1.-sin( phi(lk)*pai/180.))
*      ccr=(1.+sin(rphi(lk)*pai/180.))/(1.-sin(rphi(lk)*pai/180.))
      cc =(1.+sin( phi(lk)))/(1.-sin( phi(lk)))
      ccr=(1.+sin(rphi(lk)))/(1.-sin(rphi(lk)))


        if(ifrac(lk).eq.0.and.ps22.lt.-qu(lk)+ps11*cc)then
*        write(6,*)lk,ps1,ps2,cc,ccr
        ifrac(lk)=-2
        endif

        if(ifrac(lk).eq.-1.and.ps22.lt.-qu(lk)+ps11*cc)then
        ifrac(lk)=-4
        endif

        if(ifrac(lk).eq.-2.or.ifrac(lk).eq.-4)then
*        ps2=-rqu(lk)+ps1*ccr+alpha(lk)*elepres(lk)
      ps2=-rqu(lk)+ps1*ccr
        endif

      call sig2(ps1,ps2,0.,-alf,sx,sy,sxy)

        if(itemp.eq.1)then
        SIGM(1,LK)=sx+elepres(lk)*alpha(lk)
        SIGM(2,LK)=sy+elepres(lk)*alpha(lk)
        else
        SIGM(1,LK)=sx
        SIGM(2,LK)=sy
        endif

      SIGM(3,LK)=sxy

*      write(6,*)'ps1,ps2,pa',ps1,ps2,pa
*      write(6,*)(sigm(i,lk),i=1,3)

      DO 60 I=1,3

*      write(*,*)'sigt,sigm',i,lk,sigt(i),sigm(i,lk)

      SIGT(I)= SIGT(I)-SIGM(I,LK)
   60 continue

*      write(*,3812)(sigt(i),i=1,3)
* 3812 format(1h ,'sigt', 3g12.3)

*        if(Itemp.eq.1)then
*        PRST(1,LK)=PS1-elepres(lk)*alpha(lk)
*        PRST(2,LK)=PS2-elepres(lk)*alpha(lk)
*        else
        PRST(1,LK)=PS1
        PRST(2,LK)=PS2
*        endif

      PRST(3,LK)=alf

 3399 FORMAT(1H ,I4,2f10.3,3f10.0,5f10.2,f7.2,I5)

*        if(ifrac(lk).eq.0)goto8111

* Assign nodal load according to stress change

*        if(ifracr.eq.0.and.ifrac(lk).ne.0)then

      DO 61 I=1,6
      K=IIL(I)
      S=0.0

      DO 62 J=1,3
      S=S+BMAT(J,I,LK)*SIGT(J)
   62 continue

      S=S*VOLU(LK)
*        IF( ABS(S).GT.FOM )then
*        FOM=ABS(S)
*        KM=K
*        endif

*      write(*,*)'dump',dump

*      write(*,*)'i,k,lk,s',i,k,lk,s

      LOAD(K)=LOAD(K)+(S-sr(i,lk))*dump

        IF( ABS(S-sr(i,lk)).GT.FOM )then
        FOM=ABS(S-sr(i,lk))
        KM=K
        endif

      sr(i,lk)=s

   61 continue

*        endif

 8111 continue

*      write(6,9999)lk,(load(iil(i)),i=1,6)
* 9999 format(1h ,i4,6g12.4)

      RETURN
      END
***************************************
* PRINCIPAL STRESS
*
      subroutine sig(sx,sy,sxy,s1,s2,ph)
*      write(*,*)sx,sy,sxy
      pai=3.141593
        IF((SX-SY).EQ.0.)THEN
          IF(SXY.GT.0.)PH=PAI/4.
          IF(SXY.LT.0.)PH=-PAI/4.
          if(sxy.eq.0.)ph=0.
        ELSE
        PH=0.5*ATAN((2.*SXY)/(SX-SY))
          IF((SX-SY).LT.0.)PH=PH+PAI/2.
        ENDIF
      CP=COS(PH)
      SP=SIN(PH)
      S1=CP**2*SX+SP**2*SY+2.*SP*CP*SXY
      S2=SP**2*SX+CP**2*SY-2.*SP*CP*SXY
      PH=PH*180./PAI
      return
      END
*************************************
*
* coodinate conversion
*
      subroutine sig2(sx,sy,sxy,alph,ssx,ssy,ssxy)
      dimension s1(2,2),s2(2,2),c(2,2)
      PAI=3.141593

*      write(6,*)sx,sy,sxy,alph,ssx,ssy,ssxy

      alpha=alph*pai/180.
      sa=sin(alpha)
      ca=cos(alpha)
      s1(1,1)=sx
      s1(1,2)=sxy
      s1(2,1)=sxy
      s1(2,2)=sy
      c(1,1)= ca
      c(1,2)= sa
      c(2,1)=-sa
      c(2,2)= ca
      do 100 i=1,2
      do 100 j=1,2
      s2(i,j)=0.
  100 continue
      do 200 i=1,2
      do 200 j=1,2
      do 200 k=1,2
      do 200 l=1,2
      s2(i,l)=s2(i,l)+c(i,j)*s1(j,k)*c(l,k)
  200 continue
      ssx= s2(1,1)
      ssy= s2(2,2)
      ssxy=s2(1,2)

      return
      END
*****************************************
      SUBROUTINE ANACOL                                                 
      COMMON/C1/KKK(8),NOD(9999,3),X(9998,2),LOAD(19996),
     >          JP(999),HENI(999)    
      COMMON/C2/E(9999),V(9999),THICK(9999),EARTH(9999),VOLU(9999),FOME,     
     1          PHI(9999),qu(9999),RPHI(9999),Rqu(9999),TNS(9999),         
     2          RTNS(9999),DIL(9999),STRF(3)                              
      COMMON/C3/ST(19996,998),XMK(6,6),XX(19996),BMAT(3,6,9999),              
     1          STRESS(3,6,9999),XXT(19996)                                
      COMMON/C4/ISTAT(9999),SIGM(3,9999),STRM(3,9999),IIII,LK,FOM,KM,      
     1          ORX(9999),ORY(9999),NA,NAL
*	  COMMON/C6/a0,c0,es0,d0,s0,al0,az0,fk01,z	 
      DIMENSION SIGT(3),IIL(6)                                          
*      DOUBLE PRECISION ST
      REAL LOAD                                                         

      DO 7 I=1,3
      IIL(2*I-1)=2*NOD(LK,I)-1
      IIL(2*I)=2*NOD(LK,I)
    7 continue

      DO 17 I=1,3
      SIGT(I)=SIGM(I,LK)
   17 continue

      DO 61 I=1,6
      K=IIL(I)
      S=0.0
        DO 62 J=1,3
        S=S+BMAT(J,I,LK)*SIGT(J)
   62   continue
      S=S*VOLU(LK)
      LOAD(K)=LOAD(K)+S
   61 CONTINUE

      RETURN
      END
*
      SUBROUTINE SOLVER
      COMMON/C1/KKK(8),NOD(9999,3),X(9998,2),LOAD(19996),
     >          JP(999),HENI(999)    
      COMMON/C2/E(9999),V(9999),THICK(9999),EARTH(9999),VOLU(9999),FOME,     
     1          PHI(9999),qu(9999),RPHI(9999),Rqu(9999),TNS(9999),         
     2          RTNS(9999),DIL(9999),STRF(3)                              
      COMMON/C3/ST(19996,998),XMK(6,6),XX(19996),BMAT(3,6,9999),              
     1          STRESS(3,6,9999),XXT(19996)                                
      COMMON/C4/ISTAT(9999),SIGM(3,9999),STRM(3,9999),IIII,LK,FOM,KM,      
     1          ORX(9999),ORY(9999),NA,NAL
*	  COMMON/C6/a0,c0,es0,d0,s0,al0,az0,fk01,z	 
*      DOUBLE PRECISION SUM,TEMP,ST                                      
      REAL LOAD                                                         
C                                                                       
      NPRESC=KKK(3)                                                     
      N=KKK(5)                                                          
      JBW=KKK(6)                                                        
C   *    *    *    *    *    *    *    *    *    *    *    *    *    *  
C   *    *    *    *    *    *    *    *    *    *    *    *    *    *  
C     INTRODUCTION OF PRESCRIBED DISPLACEMENTS                          
C   *    *    *    *    *    *    *    *    *    *    *    *    *    *  
      IF( IIII ) 10,10,1000
   10 DO 100 II=1,NPRESC
      IA=JP(II)
      ST(IA,1)=ST(IA,1)*0.1E+20
*      ST(IA,1)=ST(IA,1)*0.1E+12
*      LOAD(IA)=load(ia)+ST(IA,1)*HENI(II)
      LOAD(IA)=ST(IA,1)*HENI(II)
  100 continue

*      do 9212 i=1,n/2
*      write(*,9211)i,load(i*2-1),load(i*2)
* 9212 continue
* 9211 format(1h ,"load",i5,2g12.4)

*      do 9769 i=1,n
*      write(6,9768)i,(st(i,j),j=1,8)
* 9768 format(1h ,i5,8g10.4)
* 9769 continue

* DECOMPOSITION OF BANDMATRIX

      DO 120 I=1,N
      IP=N-I+1
      IF(IP.GT.JBW) IP=JBW
      DO 120 J=1,IP
      JQ=JBW-J
      IF(JQ.GT.I-1) JQ=I-1
      SUM=ST(I,J)
      IF(JQ.EQ.0) GO TO 122
      DO 124 K=1,JQ
      IK=I-K
      JK=J+K
  124 SUM=SUM-ST(IK,K+1)*ST(IK,JK)
  122 IF(J.NE.1) GO TO 126
      IF(SUM.LE.0.0) GO TO 5000
      TEMP=1./SQRT(SUM)
      ST(I,J)=TEMP
      GO TO 120
  126 ST(I,J)=SUM*TEMP
  120 CONTINUE

*     SOLVing BANDMATRIX

 1000 continue

      do 9213 i=1,n
      xxt(i)=0.
*      xx(i)=0.
 9213 continue

      DO 140 I=1,N                                                      
      J=I-JBW+1                                                         
      IF(I+1.LE.JBW) J=1                                                
      SUM=LOAD(I)
      IF(I.EQ.1) GO TO 140                                              
      DO 142 K=J,I-1                                                    
      IA=I-K+1                                                          
  142 SUM=SUM-ST(K,IA)*XX(K)                                            
  140 XX(I)=SUM*ST(I,1)                                                 
      DO 150 I=1,N                                                      
      II=N-I+1                                                          
      J=II+JBW-1                                                        
      IF(J.GT.N) J=N                                                    
      SUM=XX(II)                                                        
      IF(I.EQ.1) GO TO 149                                              
      IB=II+1                                                           
      DO 152 K=IB,J                                                     
      IA=K-II+1                                                         
  152 SUM=SUM-ST(II,IA)*XX(K)                                           
  149 XX(II)=SUM*ST(II,1)                                               
  150 XXT(II)=XXT(II)+XX(II)                                            
    6 CONTINUE

*      do 9766 i=1,n
*      write(6,*)i,xxt(i)
* 9766 continue

      RETURN                                                            
 5000 WRITE(*,5001) I,J
 5001 FORMAT(1h ,'SUM IS NEGATIVE VALUE',2I5)
      STOP                                                              
      END                                                               
*
      SUBROUTINE FEM(XE,c1,itemp,icon,prstr1,prstr2,prstr3)
C   *    *    *    *    *    *    *    *    *    *    *    *    *    *  
C     SUBROUTINE FOR FORMATION OF ELEMENT STIFNESS AND STRESS MATRICES  
C                                                                       
C        Z = AREA OF TRIANGULAR ELEMENT                                 
C        BMAT(3,6,300) = B-MATRIX                                       
C        C1(3,3) = D-MATRIX                                             
C        XMK(6,6) = ELEMENT STIFFNESS MATRIX                            
C        STRESS(3,6) = DISPLACEMENT-STRESS MATRIX                       
C                                                                       
C   *    *    *    *    *    *    *    *    *    *    *    *    *    *  
      COMMON/C1/KKK(8),NOD(9999,3),X(9998,2),LOAD(19996),
     >          JP(999),HENI(999)    
      COMMON/C2/E(9999),V(9999),THICK(9999),EARTH(9999),VOLU(9999),FOME,     
     1          PHI(9999),qu(9999),RPHI(9999),Rqu(9999),TNS(9999),         
     2          RTNS(9999),DIL(9999),STRF(3)                              
      COMMON/C3/ST(19996,998),XMK(6,6),XX(19996),BMAT(3,6,9999),              
     1          STRESS(3,6,9999),XXT(19996)                                
      COMMON/C4/ISTAT(9999),SIGM(3,9999),STRM(3,9999),IIII,LK,FOM,KM,      
     1          ORX(9999),ORY(9999),NA,NAL                                
      COMMON/C6/a0(9999),c0(9999),es0(9999),akk(9999),bkk(9999),
     >          dkk(9999),tkk(9999),az0(9999),fk01(9999),alf
      DIMENSION ZX(3),ZY(3),XE(3,2),C1(3,3,9999),
     >          prstr1(9999),prstr2(9999),prstr3(9999),
     >          c11(3,3),c12(3,3),c13(3,3)
*      DOUBLE PRECISION ST
      REAL LOAD                                                         
C                                                                       
      INDEX=KKK(8)
      DO 400 J=1,6                                                      
      DO 400 I=1,3
      BMAT(I,J,LK)=0.0
  400 continue

      ORX(LK)=(XE(1,1)+XE(2,1)+XE(3,1))/3.0                             
      ORY(LK)=(XE(1,2)+XE(2,2)+XE(3,2))/3.0                             
      ZX(1)=XE(2,2)-XE(3,2)                                             
      ZX(2)=XE(3,2)-XE(1,2)                                             
      ZX(3)=XE(1,2)-XE(2,2)                                             
      ZY(1)=XE(3,1)-XE(2,1)                                             
      ZY(2)=XE(1,1)-XE(3,1)                                             
      ZY(3)=XE(2,1)-XE(1,1)                                             
      Z=0.5*(XE(1,1)*ZX(1)+XE(2,1)*ZX(2)+XE(3,1)*ZX(3))
	  
        if(z.le.0.)then
        write(*,*)'Element',lk,' is invalid.'
        write(6,*)'Element',lk,' is invalid.'

        endif

      ZD=2.0*Z                                                          
      BMAT(1,1,LK)=ZX(1)/ZD                                             
      BMAT(1,3,LK)=ZX(2)/ZD                                             
      BMAT(1,5,LK)=ZX(3)/ZD                                             
      BMAT(2,2,LK)=ZY(1)/ZD                                             
      BMAT(2,4,LK)=ZY(2)/ZD                                             
      BMAT(2,6,LK)=ZY(3)/ZD                                             
      BMAT(3,1,LK)=ZY(1)/ZD                                             
      BMAT(3,2,LK)=ZX(1)/ZD                                             
      BMAT(3,3,LK)=ZY(2)/ZD                                             
      BMAT(3,4,LK)=ZX(2)/ZD                                             
      BMAT(3,5,LK)=ZY(3)/ZD                                             
      BMAT(3,6,LK)=ZX(3)/ZD                                        

      IF(INDEX.EQ.0) GO TO 50                                           

C   *    *    *    *    *    *    *    *    *    *    *    *    *    *  
C     ELASTIC CONSTANTS FOR PLANE STRESS CASE                           
C   *    *    *    *    *    *    *    *    *    *    *    *    *    *  
      CC1=E(LK)/(1.-V(LK)**2)                                           
      CC2=E(LK)*V(LK)/(1.-V(LK)**2)                                     
      CC3=E(LK)*(1.-V(LK))/2./(1.-V(LK)**2)                             
      GO TO 52                                                          
C   *    *    *    *    *    *    *    *    *    *    *    *    *    *  
C     ELASTIC CONSTANTS FOR PLANE STRAIN CASE                           
C   *    *    *    *    *    *    *    *    *    *    *    *    *    *  
   50 continue
      CC1=E(LK)*(1.-V(LK))/(1.+V(LK))/(1.-2.*V(LK))                     
      CC2=E(LK)*V(LK)/(1.+V(LK))/(1.-2.*V(LK))                          
      CC3=E(LK)/2./(1.+V(LK))                                           
   52 continue
      C1(1,1,lk)=CC1                                                       
      C1(2,1,lk)=CC2                                                       
      C1(3,1,lk)=0.                                                        
      C1(1,2,lk)=CC2                                                       
      C1(2,2,lk)=CC1                                                       
      C1(3,2,lk)=0.                                                        
      C1(1,3,lk)=0.                                                        
      C1(2,3,lk)=0.                                                        
      C1(3,3,lk)=CC3                                                       

*	  write(6,*)'lk,CC1,CC2 in FEM',lk,cc1,cc2

      if (icon.eq.1)goto8001

*************************************************************
* Introduction of the constitutive equations

      pai=3.141593
*	  akk=0.008
*	  bkk=0.00005
*	  dkk=15
*	  tk=100
	  elf=(-prstr2(lk)+prstr1(lk))*dkk(lk)
      if(elf.lt.0) elf=-1*elf
	  elf1=((1.-exp(-akk(lk)*tkk(lk)))+bkk(lk)*tkk(lk))
      fk9=fk01(lk)+elf1*elf
      ad=(1.-fk9**2)*a0(lk)/pai
      af=fk9**2*a0(lk)

      do 1834 i=1,3
      do 1834 j=1,3
      c1(i,j,lk)=0.
 1834 continue

      alf=prstr3(lk)*pai/180.-pai/2.

*      alf=0.

*      write(*,*)'icon,lk,ad,af,alf',icon,lk,ad,af,alf

      C1211=ad*(atan(c0(lk)*(prstr2(lk)/es0(lk)+1.))+pai/2.)+af
      C1212=fk9*a0(lk)
      C1222=ad*(atan(c0(lk)*(prstr1(lk)/es0(lk)+1.))+pai/2.)+af

      C12(1,1)=c1211
      C12(1,2)=c1212
      C12(1,3)=0.                                                        
      C12(2,1)=c1212
      C12(2,2)=c1222
      C12(2,3)=0.                                                        
      C12(3,1)=0.                                                        
      C12(3,2)=0.
      C12(3,3)=CC3

      C11(1,1)=    cos(alf)**2
      C11(1,2)=    sin(alf)**2
      C11(1,3)=    cos(alf)*sin(alf)
      C11(2,1)=    sin(alf)**2
      C11(2,2)=    cos(alf)**2
      C11(2,3)=-   cos(alf)*sin(alf)
      C11(3,1)=-2.*cos(alf)*sin(alf)
      C11(3,2)= 2.*cos(alf)*sin(alf)
      C11(3,3)=    cos(alf)**2-sin(alf)**2

      C13(1,1)=    cos(alf)**2
      C13(1,2)=    sin(alf)**2
      C13(1,3)=-2.*cos(alf)*sin(alf)
      C13(2,1)=    sin(alf)**2
      C13(2,2)=    cos(alf)**2
      C13(2,3)= 2.*cos(alf)*sin(alf)
      C13(3,1)=    cos(alf)*sin(alf)
      C13(3,2)=-   cos(alf)*sin(alf)
      C13(3,3)=    cos(alf)**2-sin(alf)**2

      do 1833 i=1,3
      DO 1833 j=1,3
      DO 1833 K=1,3
      do 1833 m=1,3
      c1(i,m,lk)=c1(i,m,lk)+c13(i,j)*c12(j,k)*c11(k,m)
 1833 continue

 8001 continue

*********************************************************


      VOL=THICK(LK)*Z                                                   
      VOLU(LK)=VOL

      DO 510 I=1,3
      DO 511 J=1,6
      S=0.0
      DO 512 K=1,3
  512 S=S+C1(I,K,lk)*BMAT(K,J,LK)
  511 STRESS(I,J,LK)=S
  510 CONTINUE
      DO 530 I=1,6
      DO 540 J=1,6
      S=0.0
      DO 550 K=1,3
  550 S=S+BMAT(K,I,LK)*STRESS(K,J,LK)
  540 XMK(I,J)=S*VOL
  530 CONTINUE

*      IF(EARTH(LK)) 880,880,882
  882 VT=VOL*EARTH(LK)*0.33333

      DO 95 I=1,3
      JJ=NOD(LK,I)
        if(itemp.eq.1.and.icon.eq.1)then
        LOAD(2*JJ)=LOAD(2*JJ)-VT
        endif
   95 CONTINUE
  880 RETURN
      END
