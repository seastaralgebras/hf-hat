

def RasmussenParam(p,q,r,s):
    assert p >= 2*q + r, "Error: parametrization is not correct"
    assert p >= 0 and q >= 0 and r >= 0 and s >= 0, "Error: parametrization is not correct"
    assert p > s, "Error: parametrization is not correct"
    regions = []
    if p==0:
        pass
    else:
        if q==0:
            if r==0 or r==p:
                for i in range(p-1):
                    regions.append([(i+s)%p,i%p,(i+1)%p,(i+s+1)%p])
                regions.append([(s-1)%p,(p-1)%p,0%p,s%p])
            else:
                for i in range(r-1):
                    regions.append([(i+s)%p,i%p,(i+1)%p,(i+s+1)%p])
                for i in range(p-r-1):
                    regions.append([(r+i+s)%p,(r+i)%p,(r+i+1)%p,(r+i+s+1)%p])
                regions.append([(s-1)%p,(p-1)%p,0%p,s%p])
                regions.append([(r+s-1)%p,(r-1)%p,r%p,(r+1)%p])
        else:
            if r==0:
                for i in range(q-1):
                    regions.append([(2*q-1-i)%p,i%p,(i+1)%p,(2*q-2-i)%p])
                    regions.append([(i+s)%p,(2*q-i+s-1)%p,(2*q-i-2+s)%p,(i+1+s)%p])
                for i in range(p-2*q-1):
                    regions.append([(2*q+i+s)%p,(2*q+i)%p,(2*q+i+1)%p,(2*q+i+s+1)%p])
                regions.append([(s-1)%p,(p-1)%p,0%p,(2*q-1)%p,(2*q)%p,(2*q+s)%p,(2*q+s-1)%p,s%p])
                regions.append([q%p,(q-1)%p])
                regions.append([(q-1+s)%p,(q+s)%p])
            else:
                for i in range(q-1):
                    regions.append([(2*q-1-i)%p,i%p,(i+1)%p,(2*q-2-i)%p])
                    regions.append([(i+s+r)%p,(2*q-i+s+r-1)%p,(2*q-i-2+s+r)%p,(i+1+s+r)%p])
                for i in range(r-1):
                    regions.append([(2*q+i+s)%p,(2*q+i)%p,(2*q+i+1)%p,(2*q+i+s+1)%p])
                if r==p-2*q:
                        regions.append([(r+s-1)%p,(p-1)%p,0%p,(2*q-1)%p,(2*q)%p,s%p,(s-1)%p,(r+s)%p])
                else:
                    for i in range(p-2*q-r-1):
                        regions.append([(2*q+r+i+s)%p,(2*q+r+i)%p,(2*q+r+i+1)%p,(2*q+r+i+s+1)%p])
                    regions.append([0%p,(2*q-1)%p,(2*q)%p,s%p,(s-1)%p,(p-1)%p])
                    regions.append([(r+s-1)%p,(2*q+r-1)%p,(2*q+r)%p,(2*q+r+s)%p,(2*q+r+s-1)%p,(r+s)%p])
                regions.append([q%p,(q-1)%p])
                regions.append([(r+q-1+s)%p,(r+q+s)%p])
    return regions


