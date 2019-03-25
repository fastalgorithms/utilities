c
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
c 
c         This is the end of the debugging code, and the beginning
c         of the Gamma function code proper.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c
c         This file contains two user-callable subroutines: gammanew_eval,
c         and gammanew_eval_extend. The first calculates the Gamma function
c         of a real argument (either positive or negative) to full double 
c         precision. The second calculates the Gamma function of a real 
c         argument (either positive or negative) to (almost) full extended
c         precision, provided the variables are re-declared to support such
c         precision.
c
c
        subroutine gammanew_eval(x,gam)
        implicit real *8 (a-h,o-z)
        save
c
c        This subroutine evaluates the Gamma function of a real argument;
c        it produces full double precision results for both positive and 
c        negative values of the argument.
c
c                  Input parameters:
c
c  x - the argument of which the Gamma function is to be computed
c
c                  Output parameters:
c
c  gam - the Gamma function of x (how very appropriate!)
c
c
c       . . . evaluate the Gamma function if the argument is
c             on the interval [1,2]
c
        if( (x .lt. 1) .or. (x.gt. 2) ) goto 2200
c
        call gammanew_eval2(x,gam)
        return
c
 2200 continue
c
c
c       ... if the x > 2
c
        if(x .lt. 1.5) goto 3200
c
c       find x0 on the interval [1,2] with which we start,
c       and evaluate gamma at it
c
        n=x
        x0=x-n+1
        call gammanew_eval2(x0,gam)
c
c       march!
c
        do 2400 i=2,n
c
        gam=gam*(i+x0-2)
 2400 continue
c
        return
c
 3200 continue
c
c       ... if the 0 .leq. x .leq.  1
c
        if(x .lt. 0) goto 4200
c
        x0=x+1
        call gammanew_eval2(x0,gam)
        gam=gam/x
        return
c
 4200 continue
c
c       ... if x < 0
c
        x1=-x
        n=x1
        x0=x1-n
        x0=1-x0+1
c
        call gammanew_eval2(x0,gam)
c
        xx=x0-1
        do 5200 i=1,n+2
c
        gam=gam/xx
        xx=xx-1
 5200 continue
c
        return
        end        
c
c 
c 
c 
c 
        subroutine gammanew_eval_extend(x,gam)
        implicit real *16 (a-h,o-z)
        save
c
c        This subroutine evaluates the Gamma function of a real argument;
c        it produces (almost) full extended precision results for both 
c        positive and negative values of the argument, provided the 
c        variables are re-declared to support such precision.
c
c                  Input parameters:
c
c  x - the argument of which the Gamma function is to be computed
c
c                  Output parameters:
c
c  gam - the Gamma function of x (how very appropriate!)
c
c
c
c       evaluate the Gamma function if the argument is
c       on the interval [1,2]
c
        if( (x .lt. 1) .or. (x.gt. 2) ) goto 2200
c
        call gammanew_eval3(x,gam)
        return
c
 2200 continue
c
c
c       ... if the x > 2
c
        done = 1
        if(x .lt. 3*done/2) goto 3200
c
c       find x0 on the interval [1,2] with which we start,
c       and evaluate gamma at it
c
        n=x
        x0=x-n+1
        call gammanew_eval3(x0,gam)
c
c       march!
c
        do 2400 i=2,n
c
        gam=gam*(i+x0-2)
 2400 continue
c
        return
c
 3200 continue
c
c       ... if the 0 .leq. x .leq.  1
c
        if(x .lt. 0) goto 4200
c
        x0=x+1
        call gammanew_eval3(x0,gam)
        gam=gam/x
        return
c
 4200 continue
c
c       ... if x < 0
c
        x1=-x
        n=x1
        x0=x1-n
        x0=1-x0+1
c
        call gammanew_eval3(x0,gam)
c
        xx=x0-1
        do 5200 i=1,n+2
c
        gam=gam/xx
        xx=xx-1
 5200 continue
c
        return
        end        
c
c 
c 
c 
c 
        subroutine gammanew_eval3(x,gam)
        implicit real *16 (a-h,o-z)
        dimension coefs(42)
c
        data coefs/
     1  0.886226925452758013649083741670623D+00,
     2  0.161691987244425069144349421339960D-01,
     3  0.103703363422075292057509405775161D+00,
     4  -.134118505705965264609427444590074D-01,
     5  0.904033494028884643989576361780254D-02,
     6  -.242259538437044385764617259435110D-02,
     7  0.915785997143379523502915231849406D-03,
     8  -.296890121522383830093923366409561D-03,
     9  0.100928150217797671460960524521708D-03,
     *  -.336375842059855958365225717995453D-04,
     1  0.112524564378898714084403565633111D-04,
     2  -.375498601769608326351463199728310D-05,
     3  0.125284265183407932726509306293675D-05,
     4  -.417822570478480417873660426799984D-06,
     5  0.139318222887586378078540067572924D-06,
     6  -.464480612525716981485156868788579D-07,
     7  0.154844321007218425984980044056116D-07,
     8  -.516182587481302823323689251947993D-08,
     9  0.172067842623297864468229367955008D-08,
     *  -.573573438914607316487896043927398D-09,
     1  0.191193940645124209184837634463979D-09,
     2  -.637318726733661744876732845701618D-10,
     3  0.212440679157787075404568786925963D-10,
     4  -.708137779578620656802105580438353D-11,
     5  0.236046717483456631493100121338592D-11,
     6  -.786824290629355714996723594265115D-12,
     7  0.262268505662015465618725444423847D-12,
     8  -.874213989908817094617396890755866D-13,
     9  0.291499417333710813797738200069174D-13,
     *  -.971833531575606294591946449847318D-14,
     1  0.322856823583449058974923602843996D-14,
     2  -.107469659233788012193712173282348D-14,
     3  0.367876256573690494365473260918666D-15,
     4  -.123625276708537256683380049874862D-15,
     5  0.347433987117776040343629749755013D-16,
     6  -.110914923764776138469573598500169D-16,
     7  0.686389384962107642842040168449156D-17,
     8  -.245328327855757521477785550619204D-17,
     9  -.251052453492112032869491235353891D-18,
     *  0.118067019426474278738319173650695D-18,
     1  0.182746800335155441155764424855902D-18,
     2  -.642344360403561579751828672744009D-19/
c
        save
c
        n=42
        done=1
        gam=0
c
        t=(x-1-done/2)*2 
c
        tt=1
        do 1200 i=1,n
c
        gam=gam+coefs(i)*tt
        tt=tt*t
 1200 continue
c
        return
        end
c
c 
c 
c 
c 
        subroutine gammanew_eval2(x,gam) 
        implicit real *8 (a-h,o-z)
        dimension coefs(22)
c
        data coefs/
     1  0.886226925452758027D+00,0.161691987244425025D-01,
     2  0.103703363422071934D+00,-.134118505705954070D-01,
     3  0.904033494042846197D-02,-.242259538441698247D-02,
     4  0.915785994891075934D-03,-.296890120771614378D-03,
     5  0.100928168758020265D-03,-.336375903860729251D-04,
     6  0.112523678860543027D-04,-.375495650035445620D-05,
     7  0.125310464830949639D-05,-.417909902825316389D-06,
     8  0.138824572934630909D-06,-.462835109086989667D-07,
     9  0.160743202007499598D-07,-.535845567983428334D-08,
     *  0.129318662203785113D-08,-.431075842242065260D-09,
     1  0.356478271052638709D-09,-.118826785520470270D-09/
        save
c
        n=22
        done=1
        gam=0
c
        t=(x-1-done/2)*2 
c
        tt=1
        do 1200 i=1,n
c
        gam=gam+coefs(i)*tt
        tt=tt*t
 1200 continue
c
        return
        end

