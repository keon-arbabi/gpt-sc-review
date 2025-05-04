import os, re, json, base64, difflib, pymupdf, \
    traceback, concurrent.futures
from openai import OpenAI

os.chdir('projects/gpt-sc-review')

API_KEY = os.getenv('OPENAI_API_KEY')
CELL_TYPE = 'microglia'
MODEL_A = 'gpt-4.1-2025-04-14'
MODEL_B = 'gpt-4.1-mini-2025-04-14'
MAX_CTX_A = 1_047_576
MAX_CTX_B = 1_047_576
RESERVE_TOKENS_A = 32_768
RESERVE_TOKENS_B = 32_768
MAX_WORKERS = 8

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
PROMPT B (Fact-Checking & Structural Standardization)
=============================================================
SYSTEM MESSAGE:
You are a **High-Fidelity Verification and Standardization Agent**. Your
purpose is twofold:
1.  Ensure the **factual accuracy** of the provided `Initial Summary` by
    verifying every claim against the `Full PDF Text` (and associated images).
2.  Ensure the final output strictly adheres to a **Standard Output
    Structure**.

You apply corrections minimally and preserve original phrasing whenever
possible, focusing on information pertaining to the user-specified
`CELL_TYPE`.

-------------------------------------------------------------
1) INPUTS (Unchanged)
-------------------------------------------------------------
- `Full PDF Text & Images`: The complete text extracted from the
  original scientific paper and associated relevant images. This is
  the ***only*** source of truth.
- `Initial Summary`: A structured summary of the paper, focusing on
  the `CELL_TYPE`, generated by a previous process.
- `CELL_TYPE`: The specific biological cell type that is the focus
  of the summary verification.

-------------------------------------------------------------
2) CORE INSTRUCTIONS (Revised Order & Formatting Rule)
-------------------------------------------------------------

**PRIORITY 1: STRICT SOURCE GROUNDING & FACTUAL CORRECTION**
- **Verify Facts:** Systematically examine each factual statement within the
  `Initial Summary`, especially those related to the `CELL_TYPE`. Locate the
  ***exact supporting evidence*** within the `Full PDF Text` or Figures.
- **Apply Minimal Correction Principle:**
  - **Correction:** Modify statements ***only*** to align them perfectly
    with the source material. Correct inaccuracies in data, relationships,
    terminology, gene names, etc.
  - **Removal:** If a statement is entirely unsubstantiated by the source
    material, **remove it**.
  - **Preservation:** Retain the original wording and phrasing whenever
    possible *if the statement is factually correct*. Edits must be surgical.
  - **Constrained Addition:** Add information ***only*** if its absence makes
    an existing, related statement about the `CELL_TYPE` factually
    ***incomplete or misleading*** according to the source. Insert the
    ***minimum necessary*** information directly from the source.
- **No External Knowledge:** Do ***not*** use external knowledge or make
  inferences beyond what is explicitly stated in the source material.

**PRIORITY 2: STRUCTURAL STANDARDIZATION**
- **Define Standard Output Structure:** The final output *must* follow this
  exact structure and order:
    ```
    **Quick Reference** (...)
    ---
    **Detailed Summary** (...)
    [This section MUST contain the following tags/sections internally,
     in order: <metadata>...</metadata>, <methods>...</methods>,
     <findings>...</findings>, <clinical>...</clinical>]
    ---
    **Research Implications** (...)
    ```
- **Enforce Structure:** After verifying/correcting the content (Priority 1),
  examine the structure of the (potentially modified) `Initial Summary`.
  - **If** the structure **matches** the Standard Output Structure, preserve
    it exactly (including paragraph breaks, spacing, tags, etc., as per
    Priority 1).
  - **If** the structure **deviates** (e.g., missing sections, incorrect
    order, extraneous sections like preliminary `<findings>` blocks):
    - **Reorganize:** Rearrange the *verified* content blocks to fit the
      Standard Output Structure.
    - **Remove Extraneous Content:** Delete any sections or content blocks
      from the `Initial Summary` that are not part of the Standard Output
      Structure (e.g., remove preliminary `<findings>` or `<clinical>`
      blocks that appear *before* the `**Quick Reference**`).
    - **Preserve During Reorganization:** While reorganizing, preserve the
      verified phrasing, sentences, paragraph breaks within moved blocks,
      and all original XML-like tags (e.g., `<keyFinding>`,
      `<confidenceLevel>`, `<contradictionFlag>`) associated with their
      content as much as possible. Do not alter tag attributes unless
      required by factual correction (Priority 1).

-------------------------------------------------------------
3) OUTPUT SPECIFICATION (Revised)
-------------------------------------------------------------
- Generate ***only*** the full text of the revised summary.
- The output MUST adhere strictly to the **Minimal Correction Principle** for
  content and the **Standard Output Structure** for format.
- Do ***not*** include any other text, explanations, or comments about the
  changes made.

=============================================================
END OF PROMPT B
=============================================================
'''

SYSTEM_PROMPT_C = '''
=============================================================
PROMPT C (Identify Insufficient Detail)
=============================================================
SYSTEM MESSAGE:
You are a detail-oriented analysis assistant. Your task is to review a
collection of paper summaries provided in Markdown format. Each summary
is associated with a unique identifier (PID), typically found in the
header (e.g., "# summary for [PID] (...)").

Your goal is to identify summaries that offer **minimal specific detail**
about the biological `CELL_TYPE` specified by the user. "Minimal specific
detail" means the summary mentions the `CELL_TYPE` but provides very few
or no specifics about:
- Distinct subtypes or states reported in the paper.
- Key defining marker genes for those subtypes/states.
- Specific functional roles attributed to the `CELL_TYPE` or its
  subtypes in the context of the study (e.g., disease association,
  pathway involvement).
- Quantitative changes or significant findings directly related to the
  `CELL_TYPE`.

Summaries that simply state the `CELL_TYPE` was present or analysed,
without elaborating on its characteristics or findings based on the paper,
should be flagged.

-------------------------------------------------------------
INPUT:
-------------------------------------------------------------
1.  `Markdown Content`: A string containing multiple paper summaries,
    each marked with a header like "# summary for [PID] (...)".
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
END OF PROMPT C
=============================================================
'''

client = OpenAI(api_key=API_KEY)

def extract_text_and_images(pdf_path: str) -> tuple[str, list[dict]]:
    doc = pymupdf.open(pdf_path)
    full_text = []
    images_data = []
    min_dimension = 100
    max_images = 1500

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
                ext = base_image['ext']
                image_b64 = base64.b64encode(image_bytes).decode('utf-8')
                images_data.append({
                    'page': page_num + 1,
                    'b64': image_b64,
                    'ext': ext
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

def refine_summary(pdf_text: str, images_data: list[dict],
                   initial_summary: str,
                   focus_cell_type: str) -> str:

    content_list_main = []
    content_list_main.append({
        'type': 'text',
        'text': f'Full PDF Text:\n{pdf_text}'
        })
    for img_data in images_data:
        content_list_main.append({
            'type': 'text',
            'text': f'\n--- Figure Content (Page {img_data['page']}) ---'
        })
        content_list_main.append({
            'type': 'image_url',
            'image_url': {
                'url': f'data:image/{img_data['ext']};base64,{img_data['b64']}',
                'detail': 'high'
            }
        })
    content_list_main.append({
        'type': 'text',
        'text': f'\n\nInitial Summary to Fact-Check:\n{initial_summary}'
    })
    resp = client.chat.completions.create(
        model=MODEL_B,
        messages=[
            {'role': 'system', 'content': SYSTEM_PROMPT_B},
            {'role': 'user', 'content': f'CELL_TYPE: {focus_cell_type}'},
            {'role': 'user', 'content': content_list_main},
        ],
        max_tokens=RESERVE_TOKENS_B,
        temperature=0.0,
        top_p=1.0,
        frequency_penalty=0.0,
        presence_penalty=0.0,
    )
    return resp.choices[0].message.content.strip()

def calculate_similarity_score(initial: str, final: str) -> float:
    def normalize(text: str) -> list[str]:
        lines = [line.strip() for line in text.splitlines()]
        return [line for line in lines if line]
    initial_norm = '\n'.join(normalize(initial))
    final_norm = '\n'.join(normalize(final))
    matcher = difflib.SequenceMatcher(None, initial_norm, final_norm)
    return matcher.ratio()

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
                         existing_initial_summary: str | None) -> dict | None:
    try:
        print(f'[{pid}] extracting text/images...')
        pdf_text, images_data = extract_text_and_images(pdf_path)

        initial_summary = existing_initial_summary
        if initial_summary is None:
            print(f'[{pid}] generating initial summary...')
            initial_summary = generate_summary(
                pdf_text, images_data, CELL_TYPE
            )
        else:
             print(f'[{pid}] using existing initial summary.')

        if initial_summary is None:
            print(f'[{pid}] error: initial summary generation failed.')
            return None

        print(f'[{pid}] refining summary...')
        fact_checked_summary = refine_summary(
            pdf_text, images_data, initial_summary, CELL_TYPE
        )
        similarity_score = calculate_similarity_score(
            initial_summary, fact_checked_summary
        )
        print(f'[{pid}] finished (similarity: {similarity_score:.4f}).')

        return {
            'pdf': pid,
            'cell_type': CELL_TYPE,
            'initial summary': initial_summary,
            'fact-checked summary': fact_checked_summary,
            'similarity score': similarity_score
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
        if pid in summaries and summaries[pid].get('fact-checked summary'):
            print(f'[{pid}] already processed. skipping.')
            continue
        pdf_path = os.path.join(in_dir, pdf_file)
        existing_initial = summaries.get(pid, {}).get('initial summary')
        tasks_to_run.append((pid, pdf_path, existing_initial))

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
                      pids_to_exclude: set[str] | None = None):
    json_path = os.path.join(out_dir, f'{cell_type}_summaries.json')
    md_path = os.path.join(out_dir, output_filename)
    exclude_set = pids_to_exclude or set()

    print(f'\ncollecting summaries for {cell_type} into {md_path}...')
    summaries_data = load_summaries(json_path)

    if not summaries_data:
        print('no summaries found in json to collect.')
        return

    count = 0
    written_count = 0
    with open(md_path, 'w', encoding='utf-8') as md:
        for pid in sorted(summaries_data.keys()):
            data = summaries_data[pid]
            summary = data.get('fact-checked summary', '')
            if summary:
                count += 1
                if pid not in exclude_set:
                    md.write(f'# summary for {pid} ({cell_type})\n\n')
                    md.write(f'{summary}\n\n---\n\n')
                    written_count += 1

    print(f'summary collection finished for {md_path}.')
    print(f'processed {count} summaries, wrote {written_count}.')
    if exclude_set:
        print(f'excluded {len(exclude_set)} pids.')


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
                {'role': 'system', 'content': SYSTEM_PROMPT_C},
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

if __name__ == '__main__':
    all_summaries = main()

    if all_summaries:
        out_dir = 'out'
        full_md_filename = f'{CELL_TYPE}_summaries.md'
        filtered_md_filename = f'{CELL_TYPE}_summaries_filtered.md'
        full_md_path = os.path.join(out_dir, full_md_filename)

        collect_summaries(out_dir, CELL_TYPE, full_md_filename)
        insufficient_pids = identify_insufficient_summaries(
            full_md_path, CELL_TYPE)

        if insufficient_pids:
            print("\nPIDs identified with potentially insufficient detail:")
            for pid in insufficient_pids:
                print(pid)
            collect_summaries(out_dir, CELL_TYPE, filtered_md_filename,
                              pids_to_exclude=set(insufficient_pids))
        else:
            print("\nno insufficient summaries identified. "
                  "filtered file not created.")
    else:
        print("\nno summaries available to collect or filter.")