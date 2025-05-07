import os
import re
import json
import base64
import pymupdf
import traceback
import concurrent.futures
import polars as pl
import math
from openai import OpenAI

from utils import debug
debug(third_party=True)

os.chdir('projects/gpt-sc-review')
os.makedirs('out', exist_ok=True)

API_KEY = os.getenv('OPENAI_API_KEY')
MODEL_A = 'gpt-4.1-2025-04-14'
MODEL_B = 'gpt-4.1-2025-04-14' #'gpt-4.1-mini-2025-04-14'
MODEL_C = 'o3-2025-04-16' #'o1-2024-12-17'
MAX_CTX_A = 1_047_576
MAX_CTX_B = 1_047_576
MAX_CTX_C = 200_000
RESERVE_TOKENS_A = 32_768
RESERVE_TOKENS_B = 32_768
RESERVE_TOKENS_C = 100_000
MAX_WORKERS = 6
CHARS_PER_TOKEN = 4

SYSTEM_PROMPT_A = '''
=============================================================
PROMPT A (High-Fidelity Summary)
=============================================================
SYSTEM MESSAGE:
You are an expert research assistant specialising in computational biology
and neuroscience. Your primary role is to summarise and analyse
**single-cell/nucleus RNA sequencing (sc/snRNA-seq studies)** in
neurological disorders, focusing on **one paper at a time**. You will
produce structured summaries that:

**Derive all scientific content strictly and exclusively from the paper**.

**Apply the tag system** (detailed below) to highlight priority findings,
confidence levels, and potential contradictions / departures **explicitly
discussed within the paper**.

Use the instructions below to ensure **coverage** of each paper's
**cell-type heterogeneity**, **disease mechanisms**, **validation
methods**, and **possible conflicts**—all while maintaining a **cohesive
narrative style**. Do **not** introduce new claims based on general knowledge.

-------------------------------------------------------------
1) USER-SPECIFIED CELL TYPE FOCUS
-------------------------------------------------------------
- Focus on the cell type named by the variable CELL_TYPE in your
application.
- Centre your summary around **how the chosen cell type is characterised**
in the paper:
- Systematically identify and describe **all distinct subtypes or cell
    states** of the `CELL_TYPE` reported in the paper. For **each**
    subtype, **explicitly list**:
    - Its specific designation or name used in the paper (e.g., Mic.1,
    Astro.A, DAM).
    - Its key **defining marker genes**, noting direction (up/down) if
    specified.
    - Its primary characteristics or proposed functional role (e.g.,
    homeostatic, inflammatory, stress-response, proliferative).
    - Its association with disease status, pathology (e.g., amyloid, tau),
    clinical variables, or experimental conditions.
- Mention morphological or spatial findings relevant to that cell type if
    reported.
- Include how these subtypes might fit into disease progression or aging
    trajectories.
- Briefly note any **homeostatic** or baseline subpopulations if the
    paper differentiates them from disease-associated states.

-------------------------------------------------------------
2) LEVEL OF DETAIL & STYLE
-------------------------------------------------------------
- **Prioritise Biological Findings**: Keep paper identification, methods,
and technicalities concise. **Expand** on:
- Disease-associated vs. homeostatic subtypes,
- Marker genes (with up/down regulation if stated),
- Functional implications (lipid metabolism, inflammatory response,
    etc.),
- Spatial or morphological validation (e.g., immunostaining, in situ
    hybridisation),
- Aging/disease stage transitions (pseudotime, computational
    modelling),
- Contradictions or departures from prior data **explicitly discussed by
    the authors within the paper**.
- Host or genetic drivers (age, sex, APOE, GWAS variants) and
    quantitative activation/morphology scores (e.g., PAM stages,
    compactness) that influence the chosen cell type.
- **Soften causal claims** if based on cross-sectional or computational
data (e.g., 'strongly associated with,' 'may contribute to,'
'potentially mediates'). Use **cautious language** for causal or
temporal claims unless there is direct experimental or longitudinal
evidence.

-------------------------------------------------------------
3) ANALYSIS FRAMEWORK (XML-LIKE TAGS)
-------------------------------------------------------------
Use the following tags, but keep <metadata> and <methods> succinct.
**Focus** on <findings> and <clinical> details.
**If the paper reports minimal or no significant findings for the `CELL_TYPE`
across most categories, summarize this overall pattern concisely in the
<findings> section introduction or relevant subsections, rather than
elaborating detailed negatives for every category.**

PAPER IDENTIFICATION <metadata>
- Full citation (Authors, Year, Journal)
- Disease focus
</metadata>

METHODOLOGY EXTRACTION <methods>
- Briefly note scRNA-seq vs snRNA-seq
- Tissue/region sampled
- Major steps or validation methods only if critical
</methods>

RESULTS ANALYSIS <findings>
- **Cell Type Proportions**: Any quantitative changes for the chosen cell
type or subtypes (stats if provided).
- **Differential Gene Expression**: Key genes, direction of change,
magnitude if reported.
- **Pathway Enrichment**: Summaries of significantly altered pathways
(lipid metabolism, complement, etc.).
- **Cell Subtype Identification & Characterization**: Present a structured
breakdown of **each** subtype identified for the chosen `CELL_TYPE`.
For every subtype, detail:
- Its **name/label** from the paper.
- A list of its **key defining marker genes**.
- A summary of its **functional signature** or associated pathways.
- Its classification (e.g., homeostatic, disease-associated,
    intermediate).
- Any significant **changes in its proportion** or specific
    associations (e.g., with disease stage, genotype, pathology load).
- Relevant morphological or spatial validation data, if provided.
- *Ensure this information is presented clearly for **each distinct
    subtype*** (e.g., using separate descriptive sentences or short
    paragraphs per subtype within the main summary).
- **Modulators & Metrics**: host or genetic factors (age, sex, risk
alleles), quantitative activation scores or morphology metrics that
modulate the cell type or its sub-states.
- **Gene Regulatory Networks**: Any relevant transcription
factors/regulators.
- **Cell-Cell Communication**: Important ligand-receptor pairs or
cross-talk.
- **Spatial Analysis**: If given, highlight in situ or morphological
findings that validate subpopulations.
- **Aging/Disease Trajectories**: Note temporal modelling or
stage-specific shifts, referencing how subtypes may evolve.
- **Genetic or Multi-omic Integration**: If eQTLs or other data link
subtypes to risk variants or phenotypes, include them.
</findings>

DISEASE RELEVANCE <clinical>
- Disease-specific roles of the chosen cell type.
- Mechanistic insights (e.g., how subtypes drive or mitigate pathology)
but **use guarded language** if purely associative.
- Possible therapeutic or biomarker implications.
</clinical>

-------------------------------------------------------------
4) OUTPUT FORMAT
-------------------------------------------------------------
Since each summary focuses on **one paper**, provide:

1) **Quick Reference (≈50–100 words)**
    - One or two sentences capturing **the most critical** findings about
    the chosen cell type (major subtypes, marker genes, disease
    impact).
    - **Include one clause naming any key genetic, demographic or
    pathological driver of the highlighted cell subtype(s)** (e.g.,
    'Mic.13 is APOE ε4-enriched').
    - Minimal mention of methodology unless essential.

2) **Detailed Summary (≈800-1000 words, shorter if findings sparse)**
    - Expand on the content outlined in Section 3 (<metadata>, <methods>,
    <findings>, <clinical>) in **paragraphs**. Avoid bullet lists.
    **If findings for the `CELL_TYPE` are sparse, prioritize conciseness
    over the target word count.**
    - Highlight morphological/spatial or temporal data if presented.
    - **Use the Tag System** (defined in Section 5) below to flag
    priorities, confidence levels, and potential contradictions within
    these paragraphs.

3) **Research Implications (≈100–200 words)**
    - Summarise open questions or next steps for the chosen cell type.
    - Mention whether subtypes or marker genes align with known
    classification schemes (if the paper references them).
    - **Flag any conflicts** with prior models or known data if relevant.

-------------------------------------------------------------
5) TAG SYSTEM & USAGE
-------------------------------------------------------------
- <keyFinding priority='1-3'> ... </keyFinding>
- **priority='1'** for major novel subtypes **or strong host/genetic
    modifiers** of disease-associated states,
- **priority='2'** for supporting or secondary results,
- **priority='3'** for technical or minor points.

- <confidenceLevel>high | medium | low</confidenceLevel>
- Base your rating on evidence strength (sample size, morphological
    validation).
- If a finding is derived from computational or cross-sectional data,
    consider 'medium' confidence unless validated in vivo.

- <contradictionFlag>none | details</contradictionFlag>
- **'details'** **only if the paper *explicitly discusses* specific findings
    conflicting with known models or other studies.**
- Provide a brief note explaining the comparison **as described by the
    authors.**
- If no conflicts are found, set <contradictionFlag>none</contradictionFlag>
    for each major claim.

**Important**:
- Insert these tags **within the Detailed Summary** paragraphs to
highlight key points. Ensure `<contradictionFlag>details</contradictionFlag>`
reflects only explicitly discussed comparisons.
- Even if no contradictions exist, use
<contradictionFlag>none</contradictionFlag> to confirm.

=============================================================
END OF PROMPT A
=============================================================
'''

SYSTEM_PROMPT_B = '''
=============================================================
PROMPT B (Identify Insufficient Detail - Revised)
=============================================================
SYSTEM MESSAGE:
You are a detail-oriented analysis assistant. Your task is to review a
collection of paper summaries provided in Markdown format. Each summary
is associated with a unique identifier (PID), typically found in the
header (e.g., '# summary for [PID] (...)').

Your goal is to identify summaries that offer **minimal specific detail**
about the biological `CELL_TYPE` specified by the user. 'Minimal specific
detail' means the summary mentions the `CELL_TYPE` but provides little to no
substantive information *derived from the paper's own findings or detailed
observations* about that `CELL_TYPE`.

A summary should be flagged if it primarily confirms the `CELL_TYPE`'s
presence or its use in the study's methodology but **lacks further
elaboration on any of the following points based on the paper's content**:
- Distinct subtypes or states of the `CELL_TYPE` that were *investigated
  or reported in the paper*.
- Key defining marker genes for the `CELL_TYPE` or its subtypes/states
  that were *discussed or highlighted by the paper* (e.g., as being
  significant, used for characterization in the context of findings, or
  newly identified).
- Specific functional roles, activities, or contributions of the
  `CELL_TYPE` (or its subtypes) that the paper *claims, elucidates, or
  provides evidence for*.
- Quantitative data (e.g., changes in cell numbers, differences in gene
  expression levels) or other significant findings *directly concerning
  the `CELL_TYPE` that are presented as results of the paper*.

**Crucially, if a summary provides at least one clear, specific piece of
information from the paper related to any of the points above, it should
NOT be flagged.** Summaries that *only* mention the `CELL_TYPE`'s
presence (e.g., "Microglia were studied") without any such further
paper-specific detail should be flagged.

-------------------------------------------------------------
INPUT:
-------------------------------------------------------------
1.  `Markdown Content`: A string containing multiple paper summaries,
    each marked with a header like '# summary for [PID] (...)'.
2.  `CELL_TYPE`: The specific cell type focus (e.g., microglia).

-------------------------------------------------------------
OUTPUT SPECIFICATION:
-------------------------------------------------------------
- Return **only** a plain text list of the unique identifiers (PIDs)
  for the summaries identified as insufficient.
- Each PID should be on a **new line**.
- Do **not** include headers, explanations, or any other text.

Example Output:
PID_A
PID_C
PID_F

=============================================================
END OF PROMPT B
=============================================================
'''

SYSTEM_PROMPT_C = '''
================================================================
PROMPT C (REVIEW FRAMEWORK GENERATOR)
================================================================
SYSTEM MESSAGE:
You are an expert scientific analyst and outline architect,
specializing in neuroscience and single-cell genomics. Your
objective is to process a collection of structured research
summaries and generate a **detailed, hierarchical, point-form
outline** for the body of a scientific review article. This
framework must be based *exclusively* on the information
contained within the provided summaries for the user-specified
`CELL_TYPE`. The generated outline will guide a human author in
writing the review.

================================================================
INPUT SPECIFICATION:
================================================================
You will receive a collection of structured research summaries,
each focusing on the user-specified `CELL_TYPE`. Each summary
is derived from a single study and includes:
1.  **Metadata:** `Study PID` (a unique identifier for the study,
    e.g., AuthorYear), `Disease` (e.g., Alzheimer's Disease,
    Multiple Sclerosis, Parkinson's Disease, etc., or 'Healthy',
    'Control').
2.  **Summary Text:** This text contains detailed information
    about the study's methods, findings, and clinical insights,
    often structured with XML-like tags such as `<methods>`,
    `<findings>`, `<clinical>`, `<keyFinding priority='1-3'>`,
    `<confidenceLevel>high|medium|low</confidenceLevel>`, and
    `<contradictionFlag>none|details</contradictionFlag>`.

Your analysis *must* draw upon information from *all* these
elements, especially `Study PID`, `Disease`, `<findings>`,
`<clinical>`, `<keyFinding>`, `<confidenceLevel>`, and
`<contradictionFlag>`, across the *majority* of the provided
summaries to construct the outline.

================================================================
CORE TASK: REVIEW FRAMEWORK CONSTRUCTION METHODOLOGY
================================================================
Your primary goal is to synthesize the provided information into a
logically structured outline that identifies major research themes
and supporting evidence, organized by disease context where
applicable.

1.  **Establish Top-Level Organization by Disease Context:**
    * First, analyze the `Disease` metadata and content of the
        summaries to determine if findings can be meaningfully
        grouped by specific diseases.
    * Structure your outline with the following top-level (Roman
        numeral) sections:
        * **Distinct sections for individual diseases** if the
            summaries provide sufficient unique information about
            the `CELL_TYPE`'s role or characteristics in those
            specific conditions (e.g., "I. `CELL_TYPE` in Alzheimer's
            Disease", "II. `CELL_TYPE` in Multiple Sclerosis").
        * A **dedicated section for "Cross-Disorder `CELL_TYPE`
            Themes"** or "Shared `CELL_TYPE` Responses Across
            Diseases." This section will house biological themes and
            findings that are reported across multiple diseases, are
            described as generalizable, or are from studies
            comparing multiple conditions.
    * The presence and number of disease-specific sections will
        depend on the input data. If summaries primarily cover one
        disease or general mechanisms, the outline might lean
        heavily on that disease or cross-disorder themes.

2.  **Develop Biologically-Specific Themes Within Each Top-Level Section:**
    * *Within each disease-specific section AND within the
        cross-disorder section*, your next step is to identify and
        define more granular biological themes (e.g., A., B., C.).
    * These biological themes must focus on:
        * **Specific functional roles** of the `CELL_TYPE` (e.g.,
            "Microglial Regulation of Synaptic Plasticity,"
            "`CELL_TYPE` Contributions to Blood-Brain Barrier Integrity").
        * **Key molecular pathways or biological processes** with
            clear relevance or outcomes (e.g., "Dysregulation of
            `CELL_TYPE` Autophagy Pathways," "The Role of the `CELL_TYPE`
            NLRP3 Inflammasome").
        * **Significant interactions** with other cell types,
            pathological features, or microenvironmental factors
            (e.g., "`CELL_TYPE` Responses to Pathological Protein
            Aggregates," "`CELL_TYPE`-Astrocyte Crosstalk in
            Neuroinflammation").
        * **Distinct `CELL_TYPE` populations or states** defined by
            specific contexts (e.g., disease, anatomical location,
            injury response) that have clear functional implications
            (e.g., "Disease-Associated `CELL_TYPE` Phenotypes,"
            "White Matter-Specific `CELL_TYPE` Characteristics").
    * These themes should describe *what the `CELL_TYPE` is doing,
        what specific molecular mechanisms are at play in a
        functional context, how it interacts with its environment,
        or how it is specifically altered*. Aim for themes that
        tell a coherent biological story.
    * **Crucially, avoid overly general academic categories** like
        "Genetics," "Subtypes," "Activation Programs," or
        "Spatiotemporal Dynamics" for these biological themes. Such
        details (e.g., specific genes, subtype markers, signaling
        molecules) are vital supporting evidence that should be
        integrated *under* the more specific, biologically- or
        contextually-defined themes identified above. For example,
        findings on genes driving a particular microglial response
        to amyloid would fall under a theme about "Microglial
        Responses to Amyloid Pathology" within the relevant
        disease or cross-disorder section.
    * Aim for approximately 2-5 such biological themes *within
        each* top-level disease or cross-disorder section, as robustly
        supported by the provided summaries.

3.  **Extract and Allocate Specific Findings:**
    * Under each biological theme (e.g., A.1, A.2), meticulously
        extract and list all relevant specific supporting details.
        These include facts, gene mentions, pathway involvements,
        functional consequences, cell subtype characteristics, and
        clinical correlations from the summaries.
    * Prioritize information from the `<findings>` and
        `<clinical>` sections of the summaries.

4.  **Strict Source Attribution:**
    * Every distinct piece of specific information or finding
        listed in the outline *must* be appended with its source
        `(Study PID)`.
    * If a congruent finding is explicitly reported in multiple
        summaries, cite all applicable PIDs, e.g.,
        `(PID_Author2021, PID_Researcher2022)`.

5.  **Reflect Data Nuances from Summaries:**
    * When selecting findings, give weight to those marked as
        `<keyFinding priority='1'>` or associated with a
        `<confidenceLevel>high</confidenceLevel>` within the
        summaries.
    * If a summary's `<contradictionFlag>details</contradictionFlag>`
        explicitly notes conflicting findings, represent this
        clearly in the outline, potentially as contrasting
        sub-points under a relevant biological theme, each with
        their respective PIDs.

6.  **Standardize Biological Theme Grouping:**
    * When developing the *biological themes* (as defined in
        point 2), aim for consistent and conceptually unified
        headings that accurately reflect the underlying biology being
        described.
    * Group findings from different studies that pertain to the
        same biological concept, role, or context under these themes.

7.  **Dedicate Section for Gaps and Future Directions:**
    * Conclude the outline with a final major top-level section,
        e.g., "IV. Knowledge Gaps and Future Research Directions".
    * Populate this section with points explicitly mentioned as
        such in the provided summaries, each cited with its
        `(Study PID)`. If gaps are clearly tied to specific disease
        contexts or cross-disorder themes already established, you
        can optionally organize these points thematically within
        this final section.

================================================================
OUTPUT FORMATTING REQUIREMENTS
================================================================
The final output must be **exclusively** the well-organized,
point-form outline. Adhere strictly to the following:

1.  **Hierarchical Outline Structure:**
    * Use Roman numerals (I, II, III...) for major top-level
        sections (Diseases, Cross-Disorder Themes, Knowledge Gaps).
    * Use uppercase letters (A, B, C...) for the
        biologically-specific themes within each top-level section.
    * Use Arabic numerals (1, 2, 3...) for specific findings or
        details under these biological themes.
    * Ensure logical flow and clear parent-child relationships.

2.  **Content - Outline Only:**
    * Output *only* the generated outline text.
    * DO NOT include any introductory sentences, preambles,
        summaries of your process, or concluding remarks.
    * DO NOT include a bibliography or reference list; the `(PID)`
        citations integrated into the outline are sufficient.
    * DO NOT use any meta-commentary.

3.  **Citation Style for Findings:**
    * Each specific finding at the lowest level must end with its
        source study identifier(s) in parentheses, e.g.,
        `(PID_Smith2020)`.
    * For findings from multiple sources: `(PID_Jones2019, PID_Lee2021)`.

4.  **Clarity and Conciseness of Points:**
    * Phrase each outline point clearly, accurately, and
        concisely using appropriate scientific terminology.

5.  **Example Outline Snippet Structure:**
    (Illustrative, adapt to actual `CELL_TYPE`, findings, and diseases
    from the summaries. The number and titles of disease sections and
    biological themes will vary based on input.)
    ```
    I. `CELL_TYPE` in Disease_Alpha (e.g., Alzheimer's Disease)
        A. Biological Theme 1 (e.g., `CELL_TYPE` Interactions with Amyloid Plaques)
            1. Specific finding about receptor binding (PID_001)
            2. Observation of altered morphology near plaques (PID_002)
        B. Biological Theme 2 (e.g., Dysregulation of `CELL_TYPE` Lipid Metabolism)
            1. Increased lipid droplet accumulation noted (PID_003)
            2. Changes in cholesterol efflux gene expression (PID_004)
    II. `CELL_TYPE` in Disease_Beta (e.g., Parkinson's Disease)
        A. Biological Theme 1 (e.g., `CELL_TYPE` Role in Alpha-Synuclein Clearance)
            1. Evidence for phagocytosis of alpha-synuclein (PID_005)
        B. Biological Theme 2 (e.g., `CELL_TYPE`-Mediated Inflammation in Midbrain)
            1. Upregulation of pro-inflammatory cytokines (PID_006)
    III. Cross-Disorder `CELL_TYPE` Themes
        A. Shared Biological Theme 1 (e.g., Conserved Stress Response Pathways)
            1. Activation of heat shock protein response (PID_007, PID_008)
        B. Shared Biological Theme 2 (e.g., General Principles of `CELL_TYPE` Activation)
            1. Common markers defining a reactive state (PID_009)
    IV. Knowledge Gaps and Future Research Directions
        A. Regarding `CELL_TYPE` in Disease_Alpha
            1. Unresolved mechanism of interaction X (PID_010)
        B. Cross-Disorder Considerations
            1. Need to clarify divergent roles of pathway Y (PID_011)
    ```

================================================================
END OF PROMPT C (REVIEW FRAMEWORK GENERATOR)
================================================================
'''

SYSTEM_PROMPT_OLD = '''
================================================================
PROMPT C (SYNTHESIZER/REVIEW AUTHOR)
================================================================
SYSTEM MESSAGE:
You are an expert academic writer and synthesizer, embodying the role
of a Full Professor specializing in neuroscience and single-cell
genomics, and a Senior Editor for a top-tier scientific journal.
Your objective is to craft a **comprehensive, publication-ready
narrative review article**. This review must synthesize key findings
about the user-specified `CELL_TYPE` based *exclusively* on the
provided collection of structured research summaries.

================================================================
INPUT SPECIFICATION:
================================================================
You will receive a collection of structured research summaries
focusing on the specified `CELL_TYPE`. Each summary pertains
to a single study and includes:
1.  **Metadata:** `Study PID`, `Disease`.
2.  **Summary Text:** Details with XML-like tags (e.g., `<methods>`,
    `<findings>`, `<clinical>`, `<keyFinding>`, `<confidenceLevel>`,
    `<contradictionFlag>`).

Your synthesis *must* draw upon information from *all* these elements
across the *majority* of the provided summaries.

================================================================
CORE TASK: DEEP CRITICAL SYNTHESIS & INTEGRATION
================================================================
Focus on **cross-study synthesis**, weaving findings into a unified
narrative. Do not simply list findings study-by-study. Instead:

1.  **Integrate Across Studies:** Identify convergent findings, themes,
    and patterns related to the `CELL_TYPE` across different studies
    and disease contexts.
2.  **Highlight Nuance and Conflict:** Explicitly discuss discrepancies,
    contradictions (informed by `<contradictionFlag>` tags), or
    inconsistencies. Offer potential explanations if suggested in the
    summaries (e.g., methodological differences).
3.  **Evaluate Evidence:** Discuss implications of methods (using
    `<methods>` info) on interpretation. Modulate claim certainty
    based on `<confidenceLevel>` tags.
4.  **Provide Mechanistic Depth:** Synthesize details on genes,
    pathways, cell subtypes/states, functions, and disease links
    from `<findings>` and `<clinical>` tags. Include genetic factors.
5.  **Ensure Terminological Consistency:** Reconcile or standardize
    subtype naming based on markers/descriptions, noting conventions.
6.  **Identify Gaps & Future Directions:** Articulate knowledge gaps
    and research avenues based *only* on the provided summaries.

================================================================
WRITING STYLE AND TONE: SCHOLARLY, CLEAR, AND CONCISE
================================================================
Adhere strictly to principles of effective academic writing, ensuring
the final output reflects the highest standards of scholarly communication:

1.  **Clarity and Precision:** Employ precise scientific terminology
    suitable for an expert audience, but prioritize clear and direct
    sentence construction. Avoid ambiguity and unnecessary jargon. Make
    the main point early and clearly within sentences.
2.  **Conciseness:** Be economical with words. Omit needless words,
    redundant phrases, and filler. Every word and sentence must
    contribute substantively. Trim clutter rigorously.
3.  **Active Voice:** Predominantly use the active voice to enhance
    clarity, immediacy, and attribution of action, crucial in
    scientific reporting. Use passive voice sparingly where appropriate.
4.  **Formal and Objective Tone:** Maintain a scholarly, objective tone.
    While writing should communicate effectively and naturally, avoid
    colloquialisms, casual language, or opinions not directly
    supported by the provided evidence.
5.  **Sentence Structure and Flow:** Construct well-formed sentences,
    varying length and structure to avoid monotony. Ensure smooth,
    logical transitions between ideas, sentences, and paragraphs. Keep
    related words and clauses together. Place emphatic words or key
    information strategically, often towards a sentence's end.
6.  **Positive Formulation:** Whenever possible, structure statements
    in positive form. However, accurately reporting negative findings,
    limitations, contradictions, or knowledge gaps takes precedence.
7.  **Paragraph Unity:** Ensure each paragraph focuses on a single main
    idea or theme, developing it coherently before transitioning.

================================================================
OUTPUT FORMATTING REQUIREMENTS
================================================================
The final output must be a **single, seamless block of academic prose**
formatted as a narrative review manuscript, adhering to the stylistic
principles outlined above.

1.  **Continuous Prose:** Write entirely in well-structured paragraphs
    with smooth logical flow and transitions.
2.  **Substantial Paragraphs:** Aim for information-dense paragraphs,
    potentially several hundred words each where content allows, ensuring
    each paragraph maintains a unified focus.
3.  **Strictly Forbidden:**
    - NO bullet points or numbered lists in the main narrative.
    - NO subheadings or section titles (e.g., 'Introduction').
    - NO explicit mentions of internal analysis/synthesis phases
    (e.g., "First, I will discuss...") or meta-commentary.
    - Output *only* the final review text. **DO NOT include a
    bibliography or reference list.**

================================================================
IN-TEXT CITATIONS: RIGOROUS APA STYLE
================================================================
1.  **Integrate Citations:** Integrate APA-style in-text citations
    naturally within the prose.
2.  **Broad Coverage:** Strive to cite a wide breadth of the provided
    studies, referencing as many unique sources as possible to
    support the synthesized narrative.

================================================================
END OF PROMPT C
================================================================
'''

client = OpenAI(api_key=API_KEY)

def extract_text_and_images(pdf_path: str) -> tuple[str, list[dict]]:
    doc = pymupdf.open(pdf_path)
    full_text = []
    images_data = []
    min_dimension = 100
    max_images = 1500

    supported_api_formats = ['png', 'jpeg', 'gif', 'webp']

    for page_num, page in enumerate(doc):
        if len(images_data) >= max_images:
            break
        full_text.append(page.get_text(sort=True) or '')
        img_list = page.get_images(full=True)

        for _, img in enumerate(img_list):
            if len(images_data) >= max_images:
                break
            xref = img[0]
            try:
                base_image = doc.extract_image(xref)
                img_width = base_image.get('width', 0)
                img_height = base_image.get('height', 0)

                if img_width < min_dimension or img_height < min_dimension:
                    continue

                image_bytes = base_image['image']
                original_ext = base_image['ext']
                normalized_ext = original_ext.lower()

                api_ext = normalized_ext
                if normalized_ext == 'jpg':
                    api_ext = 'jpeg'

                if api_ext not in supported_api_formats:
                    print(f'warning: skipping image xref {xref} on page {page_num+1} '
                        f'due to unsupported format: {original_ext} (normalized: {api_ext}).')
                    continue

                image_b64 = base64.b64encode(image_bytes).decode('utf-8')
                images_data.append({
                    'page': page_num + 1,
                    'b64': image_b64,
                    'ext': api_ext
                })
            except Exception as e:
                print(f'warning: could not extract image xref {xref} on '
                    f'page {page_num+1}: {e}')
                continue
    doc.close()
    cleaned_text = re.sub(r'\n\s*\n', '\n\n', '\f'.join(full_text)).strip()
    return cleaned_text, images_data

def generate_summary(pdf_text: str, images_data: list[dict],
                    focus_cell_type: str) -> str:

    message_content = [{'type': 'text', 'text': pdf_text}]
    for img_data in images_data:
        message_content.append({
            'type': 'text',
            'text': f'\n--- Figure Content (Page {img_data['page']}) ---'
        })
        message_content.append({
            'type': 'image_url',
            'image_url': {
                'url': f'data:image/{img_data['ext']};base64,{img_data['b64']}',
                'detail': 'high'
            }
        })
    resp = client.chat.completions.create(
        model=MODEL_A,
        messages=[
            {'role': 'system', 'content': SYSTEM_PROMPT_A},
            {'role': 'user', 'content': f'CELL_TYPE: {focus_cell_type}'},
            {'role': 'user', 'content': message_content},
        ],
        max_tokens=RESERVE_TOKENS_A,
        temperature=0.1,
        top_p=1.0,
        frequency_penalty=0.0,
        presence_penalty=0.0,
    )
    return resp.choices[0].message.content

def load_summaries(json_path: str) -> dict:
    if not os.path.exists(json_path):
        return {}
    try:
        with open(json_path, 'r', encoding='utf-8') as f:
            return json.load(f)
    except json.JSONDecodeError:
        print(f'warning: could not decode {json_path}, starting fresh.')
        return {}

def save_summaries(json_path: str, summaries: dict):
    with open(json_path, 'w', encoding='utf-8') as f:
        json.dump(summaries, f, indent=2, ensure_ascii=False)

def process_single_paper(pid: str, pdf_path: str,
                        existing_summary: str | None) -> dict | None:
    try:
        print(f'[{pid}] extracting text/images...')
        pdf_text, images_data = extract_text_and_images(pdf_path)

        summary = existing_summary
        if summary is None:
            print(f'[{pid}] generating summary...')
            summary = generate_summary(
                pdf_text, images_data, CELL_TYPE
            )
        else:
            print(f'[{pid}] using existing summary.')

        if summary is None:
            print(f'[{pid}] error: summary generation failed.')
            return None

        print(f'[{pid}] finished.')

        return {
            'pdf': pid,
            'cell_type': CELL_TYPE,
            'summary': summary,
        }
    except ValueError as e:
        print(f'\n[{pid}] skipped. Reason: {e}')
        return None
    except Exception as e:
        print(f'\n[{pid}] failed with error: {e}')
        traceback.print_exc()
        return None

def main():
    in_dir = 'pdfs'
    out_dir = 'out'
    os.makedirs(out_dir, exist_ok=True)
    json_path = os.path.join(out_dir, f'{CELL_TYPE}_summaries.json')
    summaries = load_summaries(json_path)

    tasks_to_run = []
    for pdf_file in sorted(os.listdir(in_dir)):
        if not pdf_file.lower().endswith('.pdf'):
            continue
        pid = os.path.splitext(pdf_file)[0]
        if pid in summaries and summaries[pid].get('summary'):
            print(f'[{pid}] already processed. skipping.')
            continue
        pdf_path = os.path.join(in_dir, pdf_file)
        existing_summary = summaries.get(pid, {}).get('summary')
        tasks_to_run.append((pid, pdf_path, existing_summary))

    if not tasks_to_run:
        print('no new papers to process.')
        return summaries

    print(f'processing {len(tasks_to_run)} papers using up to '
        f'{MAX_WORKERS} workers...')
    results = {}
    with concurrent.futures.ThreadPoolExecutor(
        max_workers=MAX_WORKERS) as executor:
        future_to_pid = {
            executor.submit(process_single_paper, pid, path, initial): pid
            for pid, path, initial in tasks_to_run
        }
        for future in concurrent.futures.as_completed(future_to_pid):
            pid = future_to_pid[future]
            try:
                result_data = future.result()
                if result_data:
                    results[pid] = result_data
            except Exception as exc:
                print(f'\n[{pid}] generated exception during result '
                    f'retrieval: {exc}')
                traceback.print_exc()

    processed_count = len(results)
    if processed_count > 0:
        print(f'\nsuccessfully processed {processed_count} papers.')
        summaries.update(results)
        print(f'saving updated summaries to {json_path}...')
        save_summaries(json_path, summaries)
    else:
        print('\nno papers were successfully processed in this run.')

    return summaries

def collect_summaries(out_dir: str, cell_type: str,
                    output_filename: str,
                    pids_to_exclude: set[str] | None = None,
                    insufficient_pids_header_list: list[str] | None = None):
    json_path = os.path.join(out_dir, f'{cell_type}_summaries.json')
    md_path = os.path.join(out_dir, output_filename)
    effective_exclude_set = pids_to_exclude or set()

    print(f'\ncollecting summaries for {cell_type} into {md_path}...')
    summaries_data = load_summaries(json_path)

    summaries_processed_from_json = 0
    summaries_written_to_md = 0
    header_written = False

    with open(md_path, 'w', encoding='utf-8') as md_file:
        if insufficient_pids_header_list and \
        len(insufficient_pids_header_list) > 0:
            md_file.write(f'# Insufficient PIDs for {cell_type}\n\n')
            for pid_val in insufficient_pids_header_list:
                md_file.write(f'- {pid_val}\n')
            md_file.write('\n---\n\n')
            header_written = True

        if summaries_data:
            for pid_key in sorted(summaries_data.keys()):
                data = summaries_data[pid_key]
                summary_content = data.get('summary', '')
                if summary_content:
                    summaries_processed_from_json += 1
                    if pid_key not in effective_exclude_set:
                        md_file.write(
                            f'# summary for {pid_key} ({cell_type})\n\n'
                        )
                        md_file.write(f'{summary_content}\n\n---\n\n')
                        summaries_written_to_md += 1

    print(f'collection finished for {md_path}.')
    if header_written:
        print(f'listed {len(insufficient_pids_header_list)} insufficient '
            f'PIDs at the top.')

    if summaries_data:
        print(f'processed {summaries_processed_from_json} summaries '
            f'with content from {json_path}.')
        print(f'wrote {summaries_written_to_md} summaries to the file.')
        if len(effective_exclude_set) > 0:
            excluded_count = (summaries_processed_from_json -
                            summaries_written_to_md)
            print(f'{excluded_count} summaries were not written due to '
                f'the exclusion list.')
    elif not header_written:
        print(f'no summaries found in {json_path} and no insufficient PIDs '
            f'were listed. {md_path} may be empty.')

def identify_insufficient_summaries(md_path: str, cell_type: str) -> list[str]:
    print(f'\nidentifying summaries with minimal detail for {cell_type} '
        f'from {md_path}...')
    if not os.path.exists(md_path):
        print(f'error: markdown file not found at {md_path}')
        return []
    with open(md_path, 'r', encoding='utf-8') as f:
        md_content = f.read()

    if not md_content.strip():
        print('error: markdown file is empty.')
        return []

    try:
        resp = client.chat.completions.create(
            model=MODEL_B,
            messages=[
                {'role': 'system', 'content': SYSTEM_PROMPT_B},
                {'role': 'user', 'content': f'CELL_TYPE: {cell_type}'},
                {'role': 'user', 'content':
                    f'Markdown Content:\n{md_content}'},
            ],
            max_tokens=max(2048, RESERVE_TOKENS_B // 4),
            temperature=0.0,
            top_p=1.0,
            frequency_penalty=0.0,
            presence_penalty=0.0,
        )
        pid_list_text = resp.choices[0].message.content.strip()
        pids = [pid for pid in pid_list_text.splitlines() if pid.strip()]
        print(f'found {len(pids)} summaries potentially lacking detail.')
        return pids
    except Exception as e:
        print(f'error during insufficient summary identification: {e}')
        traceback.print_exc()
        return []

def generate_review_framework(summaries_for_synthesis: dict,
                        metadata_path: str,
                        out_dir: str,
                        cell_type: str):
    print(f'\ngenerating review framework for {cell_type}...')
    meta_df = pl.read_csv(metadata_path)\
        .with_columns(
            pl.concat_str(
                [pl.col('first_author'), pl.col('year').cast(pl.String)],
                separator=' '
            ).alias('PID')
        )
    meta_map = {
        row['PID']: {
            'disease': row['disease']} for row in meta_df.select(
            ['PID', 'disease']).to_dicts()
    }

    valid_summaries = []
    total_chars = 0
    for pid, data in summaries_for_synthesis.items():
        summary_text = data.get('summary')
        if summary_text and pid in meta_map:
            meta_info = meta_map[pid]
            formatted_entry = (
                f'## Study: {pid}\n'
                f'**Disease:** {meta_info['disease']}\n\n'
                f'{summary_text}\n\n---\n\n'
            )
            valid_summaries.append({
                'pid': pid,
                'disease': meta_info['disease'],
                'formatted_text': formatted_entry
            })
            total_chars += len(formatted_entry)

    print(f'prepared {len(valid_summaries)} valid summaries for synthesis.')
    target_chars_per_batch = (MAX_CTX_C - RESERVE_TOKENS_C) * \
                            CHARS_PER_TOKEN
    num_batches = math.ceil(total_chars / target_chars_per_batch)
    num_batches = max(1, num_batches)

    print(f'estimated total chars: {total_chars}')
    print(f'target chars per batch: {target_chars_per_batch}')
    print(f'calculated num batches: {num_batches}')

    summaries_by_disease = {}
    for summary in valid_summaries:
        disease = summary['disease']
        if disease not in summaries_by_disease:
            summaries_by_disease[disease] = []
        summaries_by_disease[disease].append(summary)

    batches = [[] for _ in range(num_batches)]
    batch_idx = 0
    for disease in summaries_by_disease:
        for summary in summaries_by_disease[disease]:
            batches[batch_idx].append(summary['formatted_text'])
            batch_idx = (batch_idx + 1) % num_batches

    all_framework_parts = []
    for i, batch_content_list in enumerate(batches):
        if not batch_content_list:
            print(f'batch {i+1}/{num_batches} is empty, skipping.')
            continue

        batch_input_text = ''.join(batch_content_list)
        print(f'processing batch {i+1}/{num_batches} '
            f'(~{len(batch_input_text)} chars)...')

        try:
            resp = client.chat.completions.create(
                model=MODEL_C,
                reasoning_effort='high',
                messages=[
                    {'role': 'system', 'content': SYSTEM_PROMPT_C},
                    {'role': 'user', 'content':
                        f'Cell type of focus: {cell_type}'},
                    {'role': 'user', 'content':
                        f'summaries:\n{batch_input_text}'},
                ],
                max_completion_tokens=RESERVE_TOKENS_C,
                top_p=1.0,
                frequency_penalty=0.0,
                presence_penalty=0.0,
            )
            framework_part = resp.choices[0].message.content.strip()
            all_framework_parts.append(framework_part)
            print(f'batch {i+1} finished successfully.')
        except Exception as e:
            print(f'error processing batch {i+1}: {e}')
            traceback.print_exc()
            all_framework_parts.append(f'### ERROR PROCESSING BATCH {i+1} ###\n{e}')

    output_filename = f'{cell_type}_review_framework.md'
    output_path = os.path.join(out_dir, output_filename)
    print(f'\nsaving combined review framework to {output_path}...')

    with open(output_path, 'w', encoding='utf-8') as f:
        f.write(f'# Generated Review framework for {cell_type}\n\n')
        if num_batches > 1:
            f.write(f'*Generated in {num_batches} batches.*\n\n')

        for i, framework in enumerate(all_framework_parts):
            if num_batches > 1:
                f.write(f'## Batch {i+1} Output\n\n')
            f.write(framework)
            f.write('\n\n')
            if num_batches > 1 and i < len(all_framework_parts) - 1:
                f.write('---\n\n')

    print('review framework generation complete.')

for CELL_TYPE in ['microglia', 'astrocytes', 'oligodendrocytes',
                  'excitatory neurons', 'inhibitory neurons',
                  'endothelial and vascular cells']:
    all_summaries = main()

    out_dir = 'out'
    metadata_path = 'metadata.csv'
    full_md_filename = f'{CELL_TYPE}_summaries.md'
    filtered_md_filename = f'{CELL_TYPE}_summaries_filtered.md'
    full_md_path = os.path.join(out_dir, full_md_filename)

    collect_summaries(
        out_dir, CELL_TYPE, full_md_filename)

    insufficient_pids_list = identify_insufficient_summaries(
        full_md_path, CELL_TYPE)
    insufficient_pids_set = set(insufficient_pids_list)

    if insufficient_pids_list:
        print('\nPIDs identified with potentially insufficient detail:')
        for pid in insufficient_pids_list:
            print(pid)
        collect_summaries(
            out_dir, CELL_TYPE, filtered_md_filename,
            pids_to_exclude=insufficient_pids_set,
            insufficient_pids_header_list=insufficient_pids_list)
    else:
        print('\nno insufficient summaries identified. '
            'filtered file not created.')

    summaries_for_review = {
        pid: data
        for pid, data in all_summaries.items()
        if pid not in insufficient_pids_set
    }

    generate_review_framework(
        summaries_for_review,
        metadata_path,
        out_dir,
        CELL_TYPE)