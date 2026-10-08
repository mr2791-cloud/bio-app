import streamlit as st

st.set_page_config(
    page_title="حاسبة النشاط الحيوي للمركبات",
    page_icon="💊",
    layout="centered"
)

st.title("💊 حاسبة التنبؤ بالنشاط الحيوي للمركبات")

st.write("أدخل خصائص المركب الكيميائي أدناه للحصول على تقييم مبدئي.")

st.markdown("---")

compound_name = st.text_input("اسم المركب:")

mw = st.number_input(
    "الوزن الجزيئي (Molecular Weight):",
    min_value=0.0,
    max_value=2000.0,
    value=300.0,
    step=1.0
)

log_p = st.number_input(
    "معامل التوزيع (LogP):",
    min_value=-10.0,
    max_value=10.0,
    value=2.5,
    step=0.1
)

if st.button("🔬 تحليل المركب", type="primary"):

    mw_valid = 150 <= mw <= 400
    logp_valid = 1.0 <= log_p <= 4.0

    st.markdown("### النتيجة:")

    if mw_valid and logp_valid:
        st.success("✅ المركب يقع ضمن النطاق المطلوب مبدئيًا.")

        st.info(
            f"اسم المركب: {compound_name if compound_name else 'غير محدد'}"
        )

    else:
        st.error("❌ المركب خارج أحد النطاقات المحددة.")

        if not mw_valid:
            st.warning(
                f"الوزن الجزيئي خارج النطاق المطلوب (150 - 400). القيمة الحالية: {mw}"
            )

        if not logp_valid:
            st.warning(
                f"LogP خارج النطاق المطلوب (1.0 - 4.0). القيمة الحالية: {log_p}"
            )
