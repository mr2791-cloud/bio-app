import streamlit as st
import requests
import pandas as pd
import plotly.graph_objects as go
from io import BytesIO

# ------------------ محاولة استيراد RDKit (اختياري) ------------------
try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors
    RDKIT_AVAILABLE = True
except Exception:
    RDKIT_AVAILABLE = False

# ------------------ إعدادات الصفحة ------------------
hide_streamlit_style = """
<style>
#MainMenu {visibility: hidden;}
footer {visibility: hidden;}
header {visibility: hidden;}
</style>
"""
st.markdown(hide_streamlit_style, unsafe_allow_html=True)

st.set_page_config(
    page_title="حاسبة النشاط الحيوي للمركبات",
    page_icon="🧪",
    layout="wide"
)

# ------------------ دوال مساعدة ------------------

PROPS = "MolecularWeight,XLogP,HBondDonorCount,HBondAcceptorCount,TPSA,RotatableBondCount,Title"

def fetch_pubchem(query, search_by="name"):
    """جلب الخصائص من PubChem."""
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/{search_by}/{query}/property/{PROPS}/JSON"
    r = requests.get(url, timeout=30)
    if r.status_code != 200:
        return None
    p = r.json()["PropertyTable"]["Properties"][0]
    return {
        "title": p.get("Title", query),
        "mw": float(p.get("MolecularWeight", 0) or 0),
        "logp": float(p.get("XLogP", 0) or 0),
        "hbd": int(p.get("HBondDonorCount", 0) or 0),
        "hba": int(p.get("HBondAcceptorCount", 0) or 0),
        "tpsa": float(p.get("TPSA", 0) or 0),
        "rot": int(p.get("RotatableBondCount", 0) or 0),
    }

def fetch_smiles(query, search_by="name"):
    """جلب SMILES من PubChem لعمل PAINS check."""
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/{search_by}/{query}/property/CanonicalSMILES/JSON"
    try:
        r = requests.get(url, timeout=30)
        if r.status_code == 200:
            return r.json()["PropertyTable"]["Properties"][0].get("CanonicalSMILES", "")
    except Exception:
        pass
    return ""

def compute_mr_from_smiles(smiles):
    """Molar Refractivity — لو RDKit متاح."""
    if RDKIT_AVAILABLE and smiles:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            return Crippen.MolMR(mol)
    return None

# ------------------ قواعد التقييم ------------------

def lipinski_rule(props):
    checks = {
        "MW ≤ 500": props["mw"] <= 500,
        "LogP ≤ 5": props["logp"] <= 5,
        "HBD ≤ 5": props["hbd"] <= 5,
        "HBA ≤ 10": props["hba"] <= 10,
    }
    violations = sum(1 for v in checks.values() if not v)
    if violations == 0:
        status, color = "مطابق تمامًا", "green"
    elif violations == 1:
        status, color = "مخالفة واحدة مسموحة", "orange"
    else:
        status, color = "مرفوض", "red"
    return checks, violations, status, color

def veber_rule(props):
    checks = {
        "Rotatable Bonds ≤ 10": props["rot"] <= 10,
        "TPSA ≤ 140 Å²": props["tpsa"] <= 140,
    }
    violations = sum(1 for v in checks.values() if not v)
    status = "مطابق" if violations == 0 else "مخالف"
    color = "green" if violations == 0 else "red"
    return checks, violations, status, color

def ghose_rule(props, mr):
    # Ghose: 160 ≤ MW ≤ 480, -0.4 ≤ LogP ≤ 5.6, 40 ≤ MR ≤ 130, atoms 20–70
    checks = {
        "160 ≤ MW ≤ 480": 160 <= props["mw"] <= 480,
        "-0.4 ≤ LogP ≤ 5.6": -0.4 <= props["logp"] <= 5.6,
        "MR بين 40 و 130": (mr is not None and 40 <= mr <= 130) if mr else None,
    }
    # عدد الذرات — تقدير من MW لو RDKit مش متاح
    atoms = None
    violations = sum(1 for v in checks.values() if v is False)
    if any(v is None for v in checks.values()):
        status, color = "غير مكتمل (RDKit غير متاح)", "gray"
    else:
        status = "مطابق" if violations == 0 else "مخالف"
        color = "green" if violations == 0 else "red"
    return checks, violations, status, color

def egan_rule(props):
    # Egan: TPSA ≤ 131.6 و -0.4 ≤ LogP ≤ 5.6
    checks = {
        "TPSA ≤ 131.6 Å²": props["tpsa"] <= 131.6,
        "-0.4 ≤ LogP ≤ 5.6": -0.4 <= props["logp"] <= 5.6,
    }
    violations = sum(1 for v in checks.values() if not v)
    status = "مطابق" if violations == 0 else "مخالف"
    color = "green" if violations == 0 else "red"
    return checks, violations, status, color

# ------------------ PAINS ------------------

# SMARTS مبسطة لمجموعات PAINS شائعة (ليست شاملة، لكنها بتغطي أشهر الأنماط)
PAINS_PATTERNS = {
    "Catechol": "c1ccccc1(O)O",
    "Quinone": "O=C1C=CC(=O)C=C1",
    "Rhodanine": "S1C(=O)NC(=S)C1",
    "Nitro-aromatic": "[N+](=O)[O-]c1ccccc1",
    "Michael acceptor (α,β-unsaturated carbonyl)": "C=CC=O",
    "Hydrazone": "C=NNC=O",
    "Azo": "N=N",
    "Isothiocyanate": "N=C=S",
}

def pains_check(smiles):
    """يرجع قائمة الأنماط المطابقة."""
    if not smiles:
        return []
    if RDKIT_AVAILABLE:
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return []
        hits = []
        for name, smarts in PAINS_PATTERNS.items():
            patt = Chem.MolFromSmarts(smarts)
            if patt and mol.HasSubstructMatch(patt):
                hits.append(name)
        return hits
    else:
        # fallback heuristic بسيط على SMILES نصيًا
        hits = []
        s = smiles.upper()
        if "C=O" in s and "C=C" in s: hits.append("Michael acceptor (تقريبي)")
        if "N=N" in s: hits.append("Azo (تقريبي)")
        if "N=C=S" in s: hits.append("Isothiocyanate (تقريبي)")
        return hits

# ------------------ Radar Chart ------------------

def make_radar(props):
    # نورم (0–1) لكل خاصية، بحيث 1 = ممتاز (داخل الحد)
    # بنستخدم نسب من الحدود الآمنة
    metrics = [
        ("MW", props["mw"], 500, "lower"),
        ("LogP", props["logp"], 5, "lower"),
        ("HBD", props["hbd"], 5, "lower"),
        ("HBA", props["hba"], 10, "lower"),
        ("TPSA", props["tpsa"], 140, "lower"),
        ("RotBonds", props["rot"], 10, "lower"),
    ]
    labels, values = [], []
    for name, val, limit, _ in metrics:
        labels.append(name)
        # قيمة معيارية: 1 لو داخل الحد، تنزل تدريجيًا لو تجاوز
        if val <= limit:
            values.append(1.0)
        else:
            values.append(max(0.0, limit / val if val else 1.0))
    # اغلاق الدائرة
    labels_closed = labels + [labels[0]]
    values_closed = values + [values[0]]

    fig = go.Figure()
    fig.add_trace(go.Scatterpolar(
        r=values_closed,
        theta=labels_closed,
        fill='toself',
        name='المركب',
        line=dict(color='royalblue'),
    ))
    # دائرة الحدود (كلها 1)
    fig.add_trace(go.Scatterpolar(
        r=[1.0] * len(labels_closed),
        theta=labels_closed,
        fill='toself',
        name='المنطقة الآمنة',
        line=dict(color='lightgreen', dash='dot'),
        opacity=0.3,
    ))
    fig.update_layout(
        polar=dict(radialaxis=dict(visible=True, range=[0, 1.2])),
        showlegend=True,
        height=500,
        margin=dict(l=40, r=40, t=40, b=40),
    )
    return fig

# ------------------ رسم وتقييم مركب ------------------

def render_compound(query, search_by, prefix=""):
    with st.spinner(f"جاري جلب بيانات {query}..."):
        props = fetch_pubchem(query, search_by)
        if not props:
            st.error(f"❌ لم يتم العثور على المركب: {query}")
            return None
        smiles = fetch_smiles(query, search_by)
        props["smiles"] = smiles
        mr = compute_mr_from_smiles(smiles)
        props["mr"] = mr

        st.success(f"✅ تم العثور على: {props['title']}")

        col_img, col_chart = st.columns([1, 1.4])
        with col_img:
            st.subheader("🖼️ التركيب البنائي")
            st.image(
                f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/{search_by}/{query}/PNG?image_size=300x300",
                caption=props["title"],
                use_column_width=True,
            )
        with col_chart:
            st.subheader("📡 Radar Chart (المنطقة الآمنة)")
            st.plotly_chart(make_radar(props), use_container_width=True)

        st.markdown("### 📊 الخصائص")
        c1, c2, c3 = st.columns(3)
        c1.metric("MW", f"{props['mw']} g/mol")
        c1.metric("LogP", props["logp"])
        c2.metric("HBD", props["hbd"])
        c2.metric("HBA", props["hba"])
        c3.metric("TPSA", f"{props['tpsa']} Å²")
        c3.metric("RotBonds", props["rot"])
        if mr:
            st.caption(f"Molar Refractivity (MR) = {mr:.2f}")

        # ----- Tabs للتقييم -----
        tab1, tab2, tab3, tab4 = st.tabs(["Lipinski", "Veber", "Ghose", "Egan"])

        with tab1:
            checks, v, status, color = lipinski_rule(props)
            for k, ok in checks.items():
                st.write(f"{'✅' if ok else '❌'} {k}")
            st.markdown(f"**النتيجة:** :{color}[{status}] (مخالفات: {v})")
            # تفسير سبب الرفض
            failed = [k for k, ok in checks.items() if not ok]
            if failed:
                st.warning(f"أسباب الرفض: {', '.join(failed)}")

        with tab2:
            checks, v, status, color = veber_rule(props)
            for k, ok in checks.items():
                st.write(f"{'✅' if ok else '❌'} {k}")
            st.markdown(f"**النتيجة:** :{color}[{status}]")

        with tab3:
            checks, v, status, color = ghose_rule(props, mr)
            for k, ok in checks.items():
                icon = "✅" if ok else ("❌" if ok is False else "⚠️")
                st.write(f"{icon} {k}")
            st.markdown(f"**النتيجة:** :{color}[{status}]")

        with tab4:
            checks, v, status, color = egan_rule(props)
            for k, ok in checks.items():
                st.write(f"{'✅' if ok else '❌'} {k}")
            st.markdown(f"**النتيجة:** :{color}[{status}]")

        # ----- PAINS -----
        st.markdown("### 🚨 PAINS Filter")
        hits = pains_check(smiles)
        if hits:
            st.error(f"⚠️ المركب قد يحتوي على أنماط PAINS: {', '.join(hits)}")
            st.caption("قد يعطي نتائج إيجابية كاذبة في الاختبارات المعملية. يُنصح بالتحقق أكثر.")
        else:
            st.success("✅ لم يتم اكتشاف أنماط PAINS شائعة.")

        return props

# ------------------ تطبيق Streamlit ------------------

st.title("🧪 حاسبة النشاط الحيوي وتصميم الأدوية")
st.write("تحليل متقدم للمركبات مع Radar Chart، فلاتر PAINS، قواعد Ghose/Egan، مقارنة، وفحص بالجملة.")

mode = st.sidebar.radio(
    "اختر الوضع:",
    ["🔬 مركب واحد", "⚖️ مقارنة مركبين", "📦 فحص بالجملة (Batch)"],
)

# =====================================================
# الوضع 1: مركب واحد
# =====================================================
if mode == "🔬 مركب واحد":
    search_type = st.radio("طريقة البحث:", ("الاسم", "SMILES"), horizontal=True)
    search_by = "name" if search_type == "الاسم" else "smiles"
    query = st.text_input("أدخل الاسم أو SMILES:", placeholder="Aspirin أو CC(=O)OC1=CC=CC=C1C(=O)O")

    if st.button("🔍 تحليل المركب", type="primary") and query.strip():
        render_compound(query.strip(), search_by)

# =====================================================
# الوضع 2: مقارنة مركبين
# =====================================================
elif mode == "⚖️ مقارنة مركبين":
    st.subheader("⚖️ مقارنة مركبين (Side-by-Side)")
    c1, c2 = st.columns(2)
    with c1:
        q1 = st.text_input("المركب الأول:", key="q1", placeholder="Aspirin")
    with c2:
        q2 = st.text_input("المركب الثاني:", key="q2", placeholder="Ibuprofen")

    if st.button("🔍 قارن", type="primary") and q1.strip() and q2.strip():
        colA, colB = st.columns(2)
        with colA:
            st.markdown(f"### 1️⃣ {q1}")
            p1 = render_compound(q1.strip(), "name")
        with colB:
            st.markdown(f"### 2️⃣ {q2}")
            p2 = render_compound(q2.strip(), "name")

        if p1 and p2:
            st.markdown("---")
            st.subheader("📋 جدول المقارنة")
            df = pd.DataFrame({
                "الخاصية": ["MW", "LogP", "HBD", "HBA", "TPSA", "RotBonds"],
                q1: [p1["mw"], p1["logp"], p1["hbd"], p1["hba"], p1["tpsa"], p1["rot"]],
                q2: [p2["mw"], p2["logp"], p2["hbd"], p2["hba"], p2["tpsa"], p2["rot"]],
            })
            st.dataframe(df, use_container_width=True)

# =====================================================
# الوضع 3: فحص بالجملة
# =====================================================
elif mode == "📦 فحص بالجملة (Batch)":
    st.subheader("📦 فحص بالجملة (حتى 50 مركب)")
    st.caption("ارفع ملف CSV/Excel فيه عمود اسمه 'name' (أو 'smiles').")
    uploaded = st.file_uploader("ارفع الملف", type=["csv", "xlsx"])

    if uploaded:
        try:
            if uploaded.name.endswith(".csv"):
                df = pd.read_csv(uploaded)
            else:
                df = pd.read_excel(uploaded)

            col = None
            for c in ["name", "Name", "compound", "smiles", "SMILES"]:
                if c in df.columns:
                    col = c
                    break
            if not col:
                st.error("❌ الملف لازم يحتوي على عمود باسم name أو smiles.")
                st.stop()

            queries = df[col].astype(str).tolist()[:50]
            search_by = "smiles" if col.lower() == "smiles" else "name"

            results = []
            progress = st.progress(0)
            for i, q in enumerate(queries):
                p = fetch_pubchem(q, search_by)
                if p:
                    smi = fetch_smiles(q, search_by)
                    hits = pains_check(smi)
                    _, lip_v, lip_s, _ = lipinski_rule(p)
                    _, veb_v, veb_s, _ = veber_rule(p)
                    results.append({
                        "Query": q,
                        "Title": p["title"],
                        "MW": p["mw"],
                        "LogP": p["logp"],
                        "HBD": p["hbd"],
                        "HBA": p["hba"],
                        "TPSA": p["tpsa"],
                        "RotBonds": p["rot"],
                        "Lipinski": lip_s,
                        "Lipinski Violations": lip_v,
                        "Veber": veb_s,
                        "Veber Violations": veb_v,
                        "PAINS": ", ".join(hits) if hits else "None",
                    })
                else:
                    results.append({"Query": q, "Title": "Not Found"})
                progress.progress((i + 1) / len(queries))

            out_df = pd.DataFrame(results)
            st.subheader("📊 نتائج الفحص")
            st.dataframe(out_df, use_container_width=True)

            # تحميل CSV
            csv_bytes = out_df.to_csv(index=False).encode("utf-8-sig")
            st.download_button(
                "⬇️ تحميل النتائج (CSV)",
                data=csv_bytes,
                file_name="batch_screening_results.csv",
                mime="text/csv",
            )
        except Exception as e:
            st.error(f"حدث خطأ: {e}")

st.markdown("---")
st.caption("⚠️ النتائج للاستخدام البحثي/التعليمي فقط، وليست بديلاً عن التجارب المعملية.")
