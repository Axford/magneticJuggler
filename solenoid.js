var DBL_EPSILON	=2.2204460492503131e-016;
var DBL_MAX 	=1.7976931348623158e+308;
var DBL_MIN 	=2.2250738585072014e-308;

var MU0 = 0.00000125663706;

function sqr(a) {
	return a * a;
}

function MAX(x,y) { return Math.max(x,y); }
function MAX3(x,y,z) { return Math.max(Math.max(x,y),z); }
function MIN(x,y) { return Math.min(x,y); }
function MIN3(x,y,z) { return Math.min(Math.min(x,y),z); }



function F(k,ierr) {
	return drf(0.0, 1.0-Math.pow(k,2), 1.0, ierr);
}


function drf( x,  y,  z, piErr)
{
    var iErr=0;
	var mu,
       xn,yn,zn,
       xndev,yndev,zndev,
       xnroot,ynroot,znroot,
       lambda,
       epsilon,
       e2,e3,
       result,
       s;
	
    var c1 = 1.0/24.0;
    var c2 = 3.0/44.0;
    var c3 = 1.0/14.0;
    var errtol = Math.pow(DBL_EPSILON*4.0,1.0/6.0);
    var lolim = 5.0*DBL_MIN;
    var hilim = DBL_MAX/5.0;

    if (piErr)
    {
        if (MIN3(x,y,z)<0.0)
        {
            iErr = 1;
        }
        else if (MIN3(x+y,x+z,y+z)<lolim)
        {
            iErr = 2;
        }
        else if (MAX3(x,y,z)>hilim)
        {
            iErr = 3;
        }
    }
    if (iErr)
    {
        if (piErr)
        {
            piErr = iErr;
        }
        result = 0.0;
    }
else
    {
        xn = x;
        yn = y;
        zn = z;

		var lc = 0;
        while (lc < 100)
        {
            mu = (xn+yn+zn)/3.0;
            xndev = 2.0-(mu+xn)/mu;
            yndev = 2.0-(mu+yn)/mu;
            zndev = 2.0-(mu+zn)/mu;
            epsilon = MAX3(Math.abs(xndev),Math.abs(yndev),Math.abs(zndev));
            if (epsilon < errtol) break;
            xnroot = Math.sqrt(xn);
            ynroot = Math.sqrt(yn);
            znroot = Math.sqrt(zn);
            lambda = xnroot*(ynroot+znroot) +ynroot*znroot;
            xn = (xn+lambda)*0.25;
            yn = (yn+lambda)*0.25;
            zn = (zn+lambda)*0.25;
            lc++;
        }
        e2 = xndev*yndev - Math.pow(zndev,2);
        e3 = xndev*yndev*zndev;
        s = 1.0 + (c1*e2-0.1-c2*e3)*e2 + c3*e3;

        if (piErr)
        {
            piErr = 0;
        }
        result = s/Math.sqrt(mu);
    }
    return result;
}


function E(k,ierr) {
	return (drf(0.0,1.0-Math.pow(k,2),1.0,ierr)
                    -(Math.pow(k,2)/3.0)*drd(0.0,1.0-Math.pow(k,2),1.0,ierr));
}

function FastE(F,k,ierr) {
	return ((F)-(Math.pow(k,2)/3.0)*drd(0.0,1.0-Math.pow(k,2),1.0,ierr));
}

function drd(x, y, z, piErr) {
    var iErr=0;
    var mu,
            xn,yn,zn,
            xndev,yndev,zndev,
            xnroot,ynroot,znroot,
            lambda,
            epsilon,
            ea,eb,ec,ed,ef,
            sigma,
            power4,
            result,
            s1,s2;
	
    var c1 = 3.0/14.0;
    var c2 = 1.0/6.0;
    var c3 = 9.0/22.0;
    var c4 = 3.0/26.0;
    var errtol = Math.pow(DBL_EPSILON/3.0,1.0/6.0);
    var uplim;
    var lolim = 2.0/Math.pow(DBL_MAX,2.0/3.0);
    var tuplim = Math.pow(DBL_MIN,1.0/3.0);
    tuplim = Math.pow(0.1*errtol,1.0/3.0)/tuplim;
    uplim = Math.pow(tuplim,2.0);

    if (piErr)
    {
        if (MIN(x,y)<0.0)
        {
            iErr = 1;
        }
        else if (MAX3(x,y,z)>uplim)
        {
            iErr = 2;
        }
        else if (MIN(x+y,z)<lolim)
        {
            iErr = 3;
        }
    }
    if (iErr)
    {
        if (piErr)
        {
            piErr = iErr;
        }
        result = 0.0;
    }
else
    {
        xn = x;
        yn = y;
        zn = z;
        sigma = 0.0;
        power4 = 1.0;
        var lc = 0;
        while (lc<100)
        {
            mu = (xn+yn+3.0*zn)*0.2;
            xndev = (mu-xn)/mu;
            yndev = (mu-yn)/mu;
            zndev = (mu-zn)/mu;
            epsilon = MAX3(Math.abs(xndev),Math.abs(yndev),Math.abs(zndev));
            if (epsilon < errtol) break;
            xnroot = Math.sqrt(xn);
            ynroot = Math.sqrt(yn);
            znroot = Math.sqrt(zn);
            lambda = xnroot*(ynroot+znroot) +ynroot*znroot;
            sigma = sigma+power4/(znroot*(zn+lambda));
            power4 = power4*0.25;
            xn = (xn+lambda)*0.25;
            yn = (yn+lambda)*0.25;
            zn = (zn+lambda)*0.25;
            lc++;
        }
        ea = xndev*yndev;
        eb = zndev*zndev;
        ec = ea-eb;
        ed = ea-6.0*eb;
        ef = ed+ec+ec;
        s1 = ed*(-c1+0.25*c3*ed-1.5*c4*zndev*ef);
        s2 = zndev*(c2*ef+zndev*(-c3*ec+zndev*c4*ea));
        if (piErr)
        {
            piErr = 0;
        }
        result = 3.0*sigma+power4*(1.0+s1+s2)/(mu*Math.sqrt(mu));
    }
    return result;
}


function testElliptics() {
	var a = 1; /* loop radius */
    var r = 0.1; /* measurement point radius */
    var x = 0; /* measurement point axial position */
    var i = 1.0; /* loop current */
    
    /* output parameters */
    var Hx,Hr; /* axial, radial field components, A/m */
    /* working vars */
    var k,q,rq,fk,ek,al,be,ga,al2,be2,alt4,Ht;
    var    ierr;

    
    /* begin computation here */
    al = r/a;
    alt4 = al*4.0;
    be = x/a;
    ga = x/r;
    al2 = al*al;
    be2 = be*be;
    q = Math.pow(1+al,2)+be2;
    rq = Math.sqrt(q);
    k = Math.sqrt(alt4/q);
    fk = F(k,ierr);
    
    ek = FastE(fk,k,ierr);
    Ht = i/(2.0*a*Math.PI*rq);
    Hx = Ht*(ek*(1-al2-be2)/(q-alt4)+fk);
    Hr = (r==0.0)?(0.0):(Ht*ga*(ek*(1+al2+be2)/(q-alt4)-fk));
    
    console.log("Axial field,  Bx: ",Hx*MU0);
    console.log("Radial field, Br: ",Hr*MU0);
}

function calcLoopField(a,r,x,i) {
	// a = loop radius
	// r = radial position
	// x = axial position
	// i = loop current

    /* output parameters */
    var Hx,Hr; /* axial, radial field components, A/m */
    /* working vars */
    var k,q,rq,fk,ek,al,be,ga,al2,be2,alt4,Ht;
    var    ierr;

    
    /* begin computation here */
    al = r/a;
    alt4 = al*4.0;
    be = x/a;
    ga = x/r;
    al2 = al*al;
    be2 = be*be;
    q = Math.pow(1+al,2)+be2;
    rq = Math.sqrt(q);
    k = Math.sqrt(alt4/q);
    fk = F(k,ierr);
    
    ek = FastE(fk,k,ierr);
    Ht = i/(2.0*a*Math.PI*rq);
    Hx = Ht*(ek*(1-al2-be2)/(q-alt4)+fk);
    Hr = (r==0.0)?(0.0):(Ht*ga*(ek*(1+al2+be2)/(q-alt4)-fk));
    
    
    // return [Axial, Radial] field components
    return [Hx,Hr];
}


function calcSolenoidField(x,y,magA,magB,magI,magLoops) {
	var B;
	var BV = new THREE.Vector2(0,0);
	
	// sum loop fields
	for (var z = 0; z < magLoops; z++) {
		var d = z * magB/magLoops - magB/2;
	
		B = calcLoopField(magA, Math.abs(x), y-d, magI);
		
		var BT = new THREE.Vector2(B[0],B[1]);
		
		BV.add(BT,BV);
		
	}
	
	return BV;
}