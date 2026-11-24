from __future__ import annotations

from fractions import Fraction


def _num(value: object, name: str) -> Fraction:
    if not isinstance(value, str) or not value or len(value)>32:
        raise ValueError(f"{name} must be a decimal string")
    try: return Fraction(value)
    except (ValueError, ZeroDivisionError) as exc: raise ValueError(f"invalid {name}") from exc


def _text(x: Fraction) -> str:
    return str(x.numerator) if x.denominator==1 else f"{x.numerator}/{x.denominator}"


def _solve3(a: list[list[Fraction]], b: list[Fraction]) -> list[Fraction]:
    matrix=[row[:] + [rhs] for row,rhs in zip(a,b)]
    for col in range(3):
        pivot=next((r for r in range(col,3) if matrix[r][col]),None)
        if pivot is None: raise ValueError("calibration system is singular")
        matrix[col],matrix[pivot]=matrix[pivot],matrix[col]
        scale=matrix[col][col]
        matrix[col]=[v/scale for v in matrix[col]]
        for row in range(3):
            if row==col: continue
            scale=matrix[row][col]
            matrix[row]=[x-scale*y for x,y in zip(matrix[row],matrix[col])]
    return [matrix[i][3] for i in range(3)]


def process(data: dict) -> dict:
    lines=data.get("reference_lines") if isinstance(data,dict) else None
    queries=data.get("query_pixels") if isinstance(data,dict) else None
    tolerance=_num(data.get("max_residual_nm"),"max_residual_nm")
    if tolerance < 0: raise ValueError("max_residual_nm must be nonnegative")
    if not isinstance(lines,list) or not 3<=len(lines)<=1000: raise ValueError("reference_lines must contain 3..1000 items")
    if not isinstance(queries,list) or not 1<=len(queries)<=10000: raise ValueError("query_pixels must contain 1..10000 items")
    ids=set(); points=[]; last_pixel=None; last_wave=None
    for line in lines:
        lid=line.get("id") if isinstance(line,dict) else None
        if not isinstance(lid,str) or not lid or lid in ids: raise ValueError("line ids must be unique")
        pixel=_num(line.get("pixel"),"pixel"); wave=_num(line.get("wavelength_nm"),"wavelength_nm")
        if pixel<0 or wave<=0: raise ValueError("pixel and wavelength ranges are invalid")
        if last_pixel is not None and (pixel<=last_pixel or wave<=last_wave): raise ValueError("reference lines must increase in both axes")
        ids.add(lid); points.append((pixel,wave,lid)); last_pixel,last_wave=pixel,wave
    sums=[sum((p[0]**k for p in points),Fraction(0)) for k in range(5)]
    rhs=[sum((p[1]*(p[0]**k) for p in points),Fraction(0)) for k in range(3)]
    coeff=_solve3([[sums[i+j] for j in range(3)] for i in range(3)],rhs)
    residuals=[]
    for pixel,wave,lid in points:
        predicted=coeff[0]+coeff[1]*pixel+coeff[2]*pixel*pixel
        residuals.append({"id":lid,"residual_nm":_text(wave-predicted)})
    worst=max((abs(Fraction(r["residual_nm"])) for r in residuals),default=Fraction(0))
    calibrated=[]
    previous=None
    for value in queries:
        pixel=_num(value,"query pixel")
        if pixel<points[0][0] or pixel>points[-1][0]: raise ValueError("query pixel is outside calibrated range")
        if previous is not None and pixel<=previous: raise ValueError("query pixels must be strictly increasing")
        previous=pixel
        wave=coeff[0]+coeff[1]*pixel+coeff[2]*pixel*pixel
        calibrated.append({"pixel":_text(pixel),"wavelength_nm":_text(wave)})
    return {"coefficients":{"c0":_text(coeff[0]),"c1":_text(coeff[1]),"c2":_text(coeff[2])},"reference_residuals":residuals,"max_abs_residual_nm":_text(worst),"status":"accepted" if worst<=tolerance else "rejected","calibrated":calibrated}

