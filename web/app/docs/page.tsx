import fs from "fs/promises";
import path from "path";
import Link from "next/link";
import ReactMarkdown from "react-markdown";

type DocEntry = {
  slug: string;
  title: string;
  content: string;
};

const formatTitle = (slug: string) =>
  slug
    .split(/[-_]/g)
    .map((chunk) => chunk.charAt(0).toUpperCase() + chunk.slice(1))
    .join(" ");

const loadDocs = async (): Promise<DocEntry[]> => {
  const candidateDirs = [
    path.join(process.cwd(), "web", "docs"),
    path.join(process.cwd(), "docs"),
  ];
  try {
    for (const dir of candidateDirs) {
      try {
        await fs.access(dir);
        const files = await fs.readdir(dir);
        const markdownFiles = files.filter((file) => file.endsWith(".md")).sort();
        if (markdownFiles.length === 0) continue;
        const docs = await Promise.all(
          markdownFiles.map(async (file) => {
            const slug = file.replace(/\.md$/, "");
            const absolutePath = path.join(dir, file);
            const raw = await fs.readFile(absolutePath, "utf8");
            const headingMatch = raw.match(/^#\s+(.*)$/m);
            const title = headingMatch?.[1]?.trim() || formatTitle(slug);
            return {
              slug,
              title,
              content: raw,
            };
          }),
        );
        return docs;
      } catch {
        continue;
      }
    }
    return [];
  } catch (error) {
    console.error("Failed to load documentation files", error);
    return [];
  }
};

const DocsPage = async ({ searchParams }: { searchParams?: { doc?: string } }) => {
  const docs = await loadDocs();
  if (docs.length === 0) {
    return (
      <div className="p-6">
        <h1 className="text-2xl font-semibold mb-4">Documentation</h1>
        <p>No documentation files were found in the <code>/docs</code> folder.</p>
      </div>
    );
  }

  const activeSlug = searchParams?.doc && docs.some((doc) => doc.slug === searchParams.doc)
    ? searchParams.doc
    : docs[0].slug;

  const activeDoc = docs.find((doc) => doc.slug === activeSlug) ?? docs[0];

  return (
    <div className="flex flex-col md:flex-row h-full w-full">
      <aside className="md:w-64 border-r border-gray-200 bg-gray-50 p-4">
        <h2 className="text-lg font-semibold mb-4">Documentation</h2>
        <nav className="space-y-2">
          {docs.map((doc) => (
            <Link
              key={doc.slug}
              href={`/docs?doc=${doc.slug}`}
              className={`block rounded px-3 py-2 text-sm ${
                doc.slug === activeSlug
                  ? "bg-black text-white"
                  : "text-gray-800 hover:bg-gray-200"
              }`}
            >
              {doc.title}
            </Link>
          ))}
        </nav>
      </aside>
      <main className="flex-1 p-6 overflow-y-auto">
        <article className="prose max-w-none">
          <ReactMarkdown>{activeDoc.content}</ReactMarkdown>
        </article>
      </main>
    </div>
  );
};

export default DocsPage;
