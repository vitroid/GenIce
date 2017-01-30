
module bond(p1=[0,0,0],p2=[1,1,1],r=1)
{
    d = p2 - p1;
    H = sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2]);
    R = sqrt(d[0]*d[0] + d[1]*d[1]);
    //atan2(y,x)
    theta = 90-atan2(d[2],R);
    phi   = atan2(d[1],d[0]);
    echo(theta,phi);
    translate(p1)
    rotate([0,theta,phi])
    cylinder(r=r, h=H);
}


bond([0,0,0], [20,40,100], r=5);

sphere(r=6);
translate([20,40,100])
    sphere(r=6);