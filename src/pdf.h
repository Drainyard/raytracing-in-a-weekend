#ifndef PDF_H
#define PDF_H

enum PDF_Type
{
    PDF_COSINE,
    PDF_HITTABLE,
    PDF_MIXTURE
};

struct PDF
{
    PDF_Type type;
    union
    {
        struct
        {
            ONB uvw;
        } cosine;
        struct
        {
            Point3 o;
            const Hittable* hittable;
        } hittable;
        struct
        {
            PDF* p[2];
        } mixture;
    };
};

PDF* alloc_pdf(PDF pdf)
{
    PDF* result = (PDF*)malloc(sizeof(PDF));
    *result = pdf;
    return result;
}

PDF cosine(Vec3& normal)
{
    PDF pdf = {};
    pdf.type = PDF_COSINE;
    pdf.cosine.uvw = build_from_w(normal);
    return pdf;
}

PDF hittable_pdf(const Hittable* hittable, const Point3& origin)
{
    PDF pdf = {};
    pdf.type = PDF_HITTABLE;
    pdf.hittable.o = origin;
    pdf.hittable.hittable = hittable;
    
    return pdf;
}

PDF mixture(PDF* p1, PDF* p2)
{
    PDF pdf = {};
    pdf.type = PDF_MIXTURE;
    pdf.mixture.p[0] = p1;
    pdf.mixture.p[1] = p2;

    return pdf;
}

f32 value(const List<Hittable>* list, const PDF& pdf, const Vec3& direction)
{
    switch(pdf.type)
    {
    case PDF_COSINE:
    {
        f32 cosine = dot(unit_vector(direction), pdf.cosine.uvw.w);
        return (cosine <= 0) ? 0 : cosine / pi;
    }
    case PDF_HITTABLE:
    {
        return pdf_value(list, *pdf.hittable.hittable, pdf.hittable.o, direction);
    }
    case PDF_MIXTURE:
    {
        return 0.5f * value(list, *pdf.mixture.p[0], direction) + 0.5f * value(list, *pdf.mixture.p[1], direction);
    }
    }
    assert(false);
    return 0.0f;
}

Vec3 generate(PDF& pdf)
{
    switch(pdf.type)
    {
    case PDF_COSINE:
    {
        return local(pdf.cosine.uvw, random_cosine_direction());
    }
    case PDF_HITTABLE:
    {
        return random(*pdf.hittable.hittable, pdf.hittable.o);
    }
    case PDF_MIXTURE:
    {
        if (random_float() < 0.5f)
        {
            return generate(*pdf.mixture.p[0]);
        }
        else
        {
            return generate(*pdf.mixture.p[1]);
        }
    }
    }
    assert(false);
    return vec3(0.0f, 0.0f, 0.0f);
}

#endif
